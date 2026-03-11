/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <boost/functional/hash.hpp>
#include <stack>
#include <utility>

#include "strpm_simd.hpp"

#define ODDFIRST 1
#define ALWAYS_RESET 0

namespace pg {

STRPM_SIMDSolver::STRPM_SIMDSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

STRPM_SIMDSolver::~STRPM_SIMDSolver()
{
}

struct RatioCompare {
    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const 
    {
        auto lhs_ratio = std::max(lhs.first, lhs.second)/std::min(lhs.first, lhs.second);
        auto rhs_ratio = std::max(rhs.first, rhs.second)/std::min(rhs.first, rhs.second);  
        
        return lhs_ratio > rhs_ratio;
    }
};

constexpr inline size_t binom(size_t n, size_t k) noexcept
{
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binom(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binom(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
}

struct ApproxSizeCompare {
    int h;

    double approxSize(int k, int t) const 
    {
        if (k == 1) return 1;

        double approximation = (k-1) * (std::log(static_cast<double>(h-1)/(k-1)) + 1) - 
                               (std::log(2* M_PI * (k-1))/2) - ((k-1)*(k-1))/static_cast<double>(2*(h-1));
    
        return (k + t)*std::log(2.0) + std::log(binom(t + k - 2, k + 2)) + approximation;
    }

    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const 
    {
        auto lhs_size = approxSize(lhs.first, lhs.second);
        auto rhs_size = approxSize(rhs.first, rhs.second);

        return lhs_size > rhs_size;
    }
};

void
STRPM_SIMDSolver::to_tmp(int idx)
{
    tmp_bits.copy_from(pm_bits[idx].data(), stdx::element_aligned);
    tmp_masks.copy_from(pm_masks[idx].data(), stdx::element_aligned);
    tmp_levels = pm_levels[idx];
}

void
STRPM_SIMDSolver::from_tmp(int idx)
{
    tmp_bits.copy_to(pm_bits[idx].data(), stdx::element_aligned);
    tmp_masks.copy_to(pm_masks[idx].data(), stdx::element_aligned);
    pm_levels[idx] = tmp_levels;
}

void
STRPM_SIMDSolver::to_best(int idx)
{
    best_bits.copy_from(pm_bits[idx].data(), stdx::element_aligned);
    best_masks.copy_from(pm_masks[idx].data(), stdx::element_aligned);
    best_levels = pm_levels[idx];
}

void
STRPM_SIMDSolver::from_best(int idx)
{
    best_bits.copy_to(pm_bits[idx].data(), stdx::element_aligned);
    best_masks.copy_to(pm_masks[idx].data(), stdx::element_aligned);
    pm_levels[idx] = best_levels;
}

void
STRPM_SIMDSolver::tmp_to_best()
{
    best_bits = tmp_bits;
    best_masks = tmp_masks;
    best_levels = tmp_levels;
}

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
STRPM_SIMDSolver::prog_tmp(int pindex, int h)
{
    // Simple case 1: Top >_p Top
    if (tmp_levels[0] == -1) return; // already Top


    simd_uint8_mask has_bits = (tmp_masks > 0);
    simd_uint8 nlb_counts { std::popcount(static_cast<uint8_t>(tmp_masks[0])) - has_bits[0] }; 
    for (size_t n = 1; n < tmp_masks.size(); n++)
    {
        // Add the NLB in this to the previous result, subtract the leading bit if there is one
        nlb_counts[n] = nlb_counts[n-1] + std::popcount(static_cast<uint8_t>(tmp_masks[n])) - has_bits[n];
    }
    
    std::vector<int>& lvls = tmp_levels;
    simd_uint8_mask smaller_than_p = simd_uint8 ([lvls](uint8_t i) { return i < lvls.size() ? lvls[i] : 0b11111111; }) <= simd_uint8{pindex};
    simd_uint8 clear_first_bit(~simd_uint8{0} & ~1);
    simd_uint8 pattern_zero_and_ones = tmp_masks & clear_first_bit;
    simd_uint8 strings_after ([lvls](uint8_t i) {
        return static_cast<uint8_t>(lvls.size() - i);
    });
    simd_uint8 needed_after ([lvls, h](uint8_t i) {
        return i < lvls.size() ? (h - lvls[i] - 1) : 0;
    });
    simd_uint8 nlb_before ([&nlb_counts](uint8_t i) { return i == 0 ? 0 : nlb_counts[i-1]; });
    simd_uint8_mask no_successor = has_bits and ((nlb_before == t) or ((nlb_counts == t) and (
             ((tmp_bits == pattern_zero_and_ones) and (strings_after == needed_after)) // Third case
          or (tmp_bits == tmp_masks))) // Fourth case
    );
    
    simd_uint8_mask has_successor = smaller_than_p and has_bits and !no_successor;

    // Use the last match to determine the successor
    int match = stdx::find_last_set(has_successor);
    if (match == -1)
    {
        if (tmp_levels[0] == 0)
        {
            // We are at top
            tmp_levels.resize(1);
            tmp_levels[0] = -1;
            return;
        }
        else 
        {
            tmp_bits[0] = 1;
            match = 0;
            tmp_levels[0] = std::min(tmp_levels[0] - 1, pindex); 
        }
    }
    else if (nlb_counts[match] == t)
    {
        // Case A: No more open bits, we erase the tail 01^j
        int reset_bits = std::countl_one(static_cast<uint8_t>(tmp_bits[match] | ~tmp_masks[match])) + 1;
        int current_bits = std::popcount(static_cast<uint8_t>(tmp_masks[match]));
        // reset the 1 to a 0
        tmp_bits[match] &= (1u << (8-reset_bits)) - 1;
        tmp_masks[match] &= (1u << (8-reset_bits)) - 1;
        // Determine new number of NLB
        // Depends on whether we reset the entire string!
        if (reset_bits < 8)
        {
            nlb_counts[match] -= current_bits - std::popcount(static_cast<uint8_t>(tmp_masks[match]));
            match++;
        }
        else 
        {
            // The new string is empty...
            nlb_counts[match] -= (current_bits - 1);
            if ((std::cmp_equal(match + 1, tmp_levels.size()) and tmp_levels[match] < h-2) or 
                (std::cmp_less(match + 1, tmp_levels.size()) and tmp_levels[match] + 1 < tmp_levels[match + 1]) )
            { 
                // Empty level!
                tmp_levels[match] = std::cmp_equal(match + 1, tmp_levels.size()) ? tmp_levels[match] + 1 : tmp_levels[match + 1] - 1;
            }
            else tmp_levels[match] ++;
        }
    }
    else 
    {
        // Case B: There are still open bits! Append to 10^j
        if ((std::cmp_less(match + 1, tmp_levels.size()) and tmp_levels[match] + 1 < tmp_levels[match + 1] and tmp_levels[match] != pindex))
        {
            // There is an empty level we can use! 
            match ++; // We append to the next string
            tmp_bits[match] = 1; 
            tmp_masks[match] = std::min(pindex, tmp_levels[match] - 1); // Move the empty string up to the next higher level
        }
        else 
        {
            int first_new = std::popcount(static_cast<uint8_t>(tmp_masks[match]));
            tmp_bits[match] |= (1u << first_new);
        }
    }
    if (match < k-1)
    {
        // Add enough bits to fill t NLB
        int bits_before = match > 0 ? nlb_counts[match - 1] : 0;
        tmp_masks[match] |= (1u << (t - bits_before + 1)) - 1;

        // Now append the new strings
        tmp_levels.resize(match+1);
        int set_to_level = tmp_levels[match] + 1;
        while (std::cmp_less(tmp_levels.size(), k-1) and set_to_level <= h-2)
        {
            tmp_levels.push_back(set_to_level);
            set_to_level++;
        }
        simd_uint8 after_match ([match](uint8_t i) { 
            return i > match; 
        });
        simd_uint8 before_levels ([lvls](uint8_t i) { 
            return i < lvls.size(); 
        });
        stdx::where((before_levels > 0 and after_match > 0), tmp_masks) = 1;
        stdx::where((after_match > 0), tmp_bits) = 0;
    }
}

/**
 * Write pm to ostream.
 */
void
STRPM_SIMDSolver::stream_pm(std::ostream &out, int idx)
{
    if (pm_levels[idx][0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int output_level = 0;
        for (int i=0; i<pm_levels.size(); i++) {
            if (i>0) out << ", ";
            while (pm_levels[idx][i] > output_level)
            {
                out << "ε, ";
                output_level++;
            }
            out << std::bitset<8>(static_cast<uint8_t>(pm_bits[idx][i])) << "/" << std::bitset<8>(static_cast<uint8_t>(pm_masks[idx][i]));
        }
        while (output_level < h-2)
        {
            out << ", ε";
            output_level++;
        }
        out << " }";
    }
}

/**
 * Write tmp to ostream.
 */
void
STRPM_SIMDSolver::stream_simd(std::ostream &out, simd_uint8& bits, simd_uint8& masks, std::vector<int>& levels)
{
    if (levels[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int output_level = 0;
        for (int i=0; i<levels.size(); i++) {
            if (i>0) out << ", ";
            while (levels[i] > output_level)
            {
                out << "ε, ";
                output_level++;
            }
            out << std::bitset<8>(static_cast<uint8_t>(bits[i])) << "/" << std::bitset<8>(static_cast<uint8_t>(masks[i]));
        }
        while (output_level < h-2)
        {
            out << ", ε";
            output_level++;
        }
        out << " }";
    }
}

/**
 * Compare tmp and best
 * res := -1 :: tmp < best
 * res := 0  :: tmp = best
 * res := 1  :: tmp > best
 */
int
STRPM_SIMDSolver::compare(int pindex)
{
    // cases involving Top
    if (tmp_levels[0] == -1 and best_levels[0] == -1) return 0;
    if (tmp_levels[0] == -1) return 1;
    if (best_levels[0] == -1) return -1;

    simd_uint8 shorter_string = tmp_masks & best_masks;
    simd_uint8 first_length_difference = (tmp_masks ^ best_masks) & ~((tmp_masks ^ best_masks) * 2);
    simd_uint8 bit_xor = tmp_bits ^ best_bits;
    simd_uint8_mask a_less = ((tmp_masks < best_masks) and (bit_xor & (shorter_string + first_length_difference)) == first_length_difference) or
                             ((tmp_masks > best_masks) and (bit_xor & (shorter_string + first_length_difference)) == 0);
    simd_uint8_mask b_less = ((best_masks < tmp_masks) and (bit_xor & (shorter_string + first_length_difference)) == first_length_difference) or
                             ((best_masks > tmp_masks) and (bit_xor & (shorter_string + first_length_difference)) == 0);


    simd_uint8 different_bits = (shorter_string & bit_xor);
    simd_uint8 first_bit_difference ([&different_bits](uint8_t i){ 
        return 1u << (std::countr_zero(static_cast<uint8_t>(different_bits[i])));
    });
    simd_uint8_mask a_greater = (different_bits > 0) and ((tmp_bits & first_bit_difference) > 0);
    simd_uint8_mask b_greater = (different_bits > 0) and ((best_bits & first_bit_difference) > 0);

    for (size_t i = 0; i < std::max(tmp_levels.size(), best_levels.size()); i++)
    {
        // One of the two has more nonempty strings but we were equal up until here
        if (i >= best_levels.size())
        {
            return (tmp_bits[i] & 1) == 0 ? -1 : 1;
        }
        else if (i >= tmp_levels.size())
        {
            return (best_bits[i] & 1) == 0 ? 1 : -1;
        }
        // We are past the considered index!
        else if (tmp_levels[i] > pindex and best_levels[i] > pindex)
        {
            break;
        }
        // One of the two has a bitstring "earlier"
        else if ((tmp_levels[i] <= pindex and best_levels[i] > pindex) or (tmp_levels[i] < best_levels[i]))
        {
            return (tmp_bits[i] & 1) == 0 ? -1 : 1;
        }
        else if ((tmp_levels[i] > pindex and best_levels[i] <= pindex) or (tmp_levels[i] > best_levels[i]))
        {
            return (best_bits[i] & 1) == 0 ? 1 : -1;
        }
        // The two levels are equal... We have to compare strings
        else if (a_greater[i] or (!a_less[i] and b_less[i]))
        {
            return 1;
        }
        else if (b_greater[i] or (a_less[i] and !b_less[i]))
        {
            return -1;
        }
    }
    
    // The birstrings are entirely equal
    return 0;
}

bool
STRPM_SIMDSolver::lift(int v, int target, int &str, int pl)
{
    // check if already Top
    if (pm_levels[v][0] == -1) return false; // already Top

    const int pr = priority(v);
    const int pindex = pl == 0 ? (h-1)-(pr+1)/2-1 : (h-1)-pr/2-1;

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[37;1mupdating vertex " << label_vertex(v) << " (" << pr << " => " << pindex << ")" << (owner(v)?" (odd)":" (even)") << "\033[m with current measure";
        stream_pm(logger, v);
        logger << std::endl;
    }
#endif

    // if even owns and target is set, just check if specific target is better
    if (owner(v) == pl and target != -1) {
        to_tmp(target);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "to target " << label_vertex(target) << "(" << target << ")" << ":";
                stream_simd (logger, tmp_bits, tmp_masks, tmp_levels);
                logger << " =>";
            }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        //else trunc_tmp(pindex);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_simd (logger, tmp_bits, tmp_masks, tmp_levels);
                logger << std::endl;
            }
#endif
        to_best(v);
        if (compare(pindex) > 0) {
            from_tmp(v);
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[32;1mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
                stream_simd (logger, tmp_bits, tmp_masks, tmp_levels);
                logger << " (to " << label_vertex(target) << ")\n";
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // compute best measure
    bool first = true;
    for (auto curedge = outs(v); *curedge != -1; curedge++) {
        int to = *curedge;
        if (disabled[to]) continue;
        to_tmp(to);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "to successor " << label_vertex(to) << " from";
            stream_simd (logger, tmp_bits, tmp_masks, tmp_levels);
            logger << " =>";
        }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        //else trunc_tmp(pindex);
#ifndef NDEBUG
        if (trace >= 2) {
            stream_simd (logger, tmp_bits, tmp_masks, tmp_levels);
            logger << std::endl;
        }
#endif
        if (first) {
            tmp_to_best();
            str = to;
        } else if (owner(v) == pl) {
            // we want the max!
            if (compare(pindex) > 0) {
                tmp_to_best();
                str = to;
            }
        } else {
            // we want the min!
            if (compare(pindex) < 0) {
                tmp_to_best();
                str = to;
            }
        }
        first = false;
    }

    // set best to pm if higher
    to_tmp(v);
    if (compare(pindex) < 0) {
#ifndef NDEBUG
        if (trace >= 1) {
            logger << "\033[1;32mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
            stream_simd (logger, best_bits, best_masks, best_levels);
            logger << " (to " << label_vertex(str) << ")\n";
        }
#endif
        from_best(v);
        return true;
    } else {
        return false;
    }
}

static int
ceil_log2(unsigned long long x)
{
    static const unsigned long long t[6] = {
        0xFFFFFFFF00000000ull,
        0x00000000FFFF0000ull,
        0x000000000000FF00ull,
        0x00000000000000F0ull,
        0x000000000000000Cull,
        0x0000000000000002ull
    };

    int y = (((x & (x - 1)) == 0) ? 0 : 1);
    int j = 32;
    int i;

    for (i = 0; i < 6; i++) {
        int k = (((x & t[i]) == 0) ? 0 : j);
        y += k;
        x >>= k;
        j >>= 1;
    }

    return y;
}

static int
floor_log2 (unsigned long long x)
{
    static const unsigned long long t[6] = {
        0xFFFFFFFF00000000ull,
        0x00000000FFFF0000ull,
        0x000000000000FF00ull,
        0x00000000000000F0ull,
        0x000000000000000Cull,
        0x0000000000000002ull
    };

    int y = 0;             // no +1 for non-powers of two[2][1]
    int j = 32;

    for (int i = 0; i < 6; i++) {
        int k = (((x & t[i]) == 0) ? 0 : j);
        y += k;
        x >>= k;
        j >>= 1;
    }

    return y;
}

struct Node
{
    int k;
    int t;
    int h;
    bool isU;
};


// Keep track of already computed sizes, this is cached beyond single calls of
// tree_size below
std::unordered_map<std::tuple<int, int, int>, unsigned,
                   boost::hash<std::tuple<int, int, int>>> treeU_simd;
std::unordered_map<std::tuple<int, int, int>, unsigned,
                   boost::hash<std::tuple<int, int, int>>> treeV_simd;

unsigned tree_size_simd(int k, int t, int h) 
{
    std::stack<Node> stack;

    stack.push ({k, t, h, true});
    while (!stack.empty())
    {
        Node& tos = stack.top();
        if (tos.isU and tos.h == 1 and tos.k == 1)
        {
            treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = 1;
            stack.pop ();
        }
        else if (tos.isU and tos.h > 1 and tos.k == 1)
        {
            auto son = treeU_simd.find(std::make_tuple(tos.k, tos.t, tos.h - 1)); 
            if (son != treeU_simd.end())
            {
                treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k, tos.t, tos.h - 1, true});
        }
        else if (tos.h >= tos.k and tos.k >= 2 and tos.t == 0)
        {
            auto son = treeU_simd.find(std::make_tuple(tos.k - 1, tos.t, tos.h - 1));
            if (son != treeU_simd.end())
            {
                if (tos.isU) treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                else treeV_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k - 1, tos.t, tos.h - 1, true});
        }
        else if (!tos.isU and tos.h >= tos.k and tos.k >= 2 and tos.t >= 1)
        {
            auto son1 = treeV_simd.find(std::make_tuple(tos.k, tos.t - 1, tos. h));
            auto son2 = treeU_simd.find(std::make_tuple(tos.k - 1, tos.t, tos.h - 1));
            if (son1 != treeV_simd.end() and son2 != treeU_simd.end())
            {
                treeV_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son1->second * 2 + son2->second;
                stack.pop ();
            }
            else
            {
                stack.push ({tos.k - 1, tos.t, tos.h - 1, true});
                stack.push ({tos.k, tos.t - 1, tos.h, false});
            }
        }
        else if (tos.isU and tos.h == tos.k and tos.k >= 2)
        {
            auto son = treeV_simd.find(std::make_tuple(tos.k, tos.t, tos.h));
            if (son != treeV_simd.end()) 
            {
                treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son->second;
                stack.pop ();
            }
            else stack.push ({tos.k, tos.t, tos.h, false});
        }
        else if (tos.isU and tos.h > tos.k and tos.k >= 2)
        {
            auto son1 = treeV_simd.find(std::make_tuple(tos.k, tos.t, tos.h));
            auto son2 = treeU_simd.find(std::make_tuple(tos.k, tos.t, tos.h - 1));
            if (son1 != treeV_simd.end() and son2 != treeU_simd.end())
            {
                treeU_simd[std::make_tuple(tos.k, tos.t, tos.h)] = son1->second * 2 + son2->second;
                stack.pop ();
            }
            else
            {
                stack.push ({tos.k, tos.t, tos.h - 1, true});
                stack.push ({tos.k, tos.t, tos.h, false});
            }
        }
        else assert(false); // We should never get here
    }
    
    return treeU_simd[std::make_tuple(k, t, h)];
}

struct SizeCompare
{
    int h;

    bool operator()(const std::pair<int,int>& lhs,
                    const std::pair<int,int>& rhs) const 
    {
        auto lhs_size = tree_size_simd(lhs.first, lhs.second, h);
        auto rhs_size = tree_size_simd(rhs.first, rhs.second, h);
        
        return lhs_size > rhs_size;
    }
};


void
STRPM_SIMDSolver::run(int t_val, int k_val, int depth, int player)
{
    // Marcin's word: think of h as the number of priorities of the
    // opponent... PLUS ONE!
    t = t_val;
    h = depth + 1;  // FIXME: This is Guillermo's hack, the +1
    k = k_val;  // Maybe possible: std::min(t + 2, h);

    tmp_levels.reserve(k-1);
    best_levels.reserve(k-1);

#ifndef NDEBUG
    logger << "Strahler-tree parameters for player " << player << ": k = " << k << ", t = " << t << ", h = " << h << std::endl;
#endif

    // initialize progress measures - Every node is set to the smallest leaf in the tree
    pm_bits = std::vector<std::vector<uint8_t>> (nodecount(), std::vector<uint8_t>(8, 0));
    std::vector<uint8_t> initial_mask (8, 0);
    initial_mask[0] = (1 << (t+1)) - 1;
    std::vector<int> initial_levels (k-1);
    for (size_t i = 1; i < k-1; i++)
    {
        initial_levels[i] = i;
        initial_mask[i] = 1;
    }
    pm_masks = std::vector<std::vector<uint8_t>> (nodecount(), initial_mask);
    pm_levels = std::vector<std::vector<int>> (nodecount(), initial_levels);

#ifndef NDEBUG
    if (trace >= 1)
    {
        logger << "Initial PM: " << std::endl;
        stream_pm(logger, 0);
        logger << std::endl;
    }
#endif

    for (int n=nodecount()-1; n>=0; n--) {
        if (disabled[n]) continue;
        lift_attempt++;
        int s;
        if (lift(n, -1, s, player)) {
            lift_count++;
            // lift_counters[n]++;
            for (auto curedge = ins(n); *curedge != -1; curedge++) {
                int from = *curedge;
                if (disabled[from]) continue;
                lift_attempt++;
                int s;
                if (lift(from, n, s, player)) {
                    lift_count++;
                    // lift_counters[from]++;
                    todo_push(from);
                }
            }
        }
    }

    while (!Q.empty()) {
        int n = todo_pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            if (disabled[from]) continue;
            lift_attempt++;
            int s;
            if (lift(from, n, s, player)) {
                lift_count++;
                // lift_counters[from]++;
                todo_push(from);
            }
        }
    }

    /**
     * Derive strategies.
     */

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm_levels[v][0] != -1) {
            if (owner(v) != player) {
                // TODO: don't rely on the strategy array in the Game class
                if (lift(v, -1, game.getStrategy()[v], player)) logger << "error: " << v << " is not progressive!" << std::endl;
            }
        }
    }

    if (trace) {
        for (int v=0; v<nodecount(); v++) {
            if (disabled[v]) continue;

            logger << "\033[1m" << label_vertex(v) << (owner(v)?" (odd)":" (even)") << "\033[m:";
            stream_pm(logger, v);

            if (pm_levels[v][0] != -1) {
                if (owner(v) != player) {
                    logger << " => " << label_vertex(game.getStrategy(v));
                }
            }

            logger << std::endl;
        }
    }

    /**
     * Mark solved.
     */

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm_levels[v][0] != -1) Solver::solve(v, 1-player, game.getStrategy(v));
    }

    Solver::flush();
}

void
STRPM_SIMDSolver::run()
{
    int max_prio = priority(nodecount()-1);

    // compute ml (max l) and the h for even/odd
    int t_max = floor_log2(nodecount());
    int h0 = (max_prio/2)+1;
    int h1 = (max_prio+1)/2;

    int h_max = std::max(h0, h1);
    int k_max = std::min(t_max + 2, h_max);

    // create datastructures
    Q.resize(nodecount());
    dirty.resize(nodecount());

    // Create a priority queue for (k, t) pairs and push init with (1, 1)
    std::priority_queue<
        std::pair<int,int>,
        std::vector<std::pair<int,int>>,
        RatioCompare
        //ApproxSizeCompare
    > pq { };
    pq.push({2, 1});
    /*
    To use SizeCompare:
    std::priority_queue<
        std::pair<int,int>,
        std::vector<std::pair<int,int>>,
        SizeCompare
    > pq { SizeCompare { h_max } };
    */

    // Keep track of already tried combinations
    std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> already_tried;

#ifndef NDEBUG
    logger << "Max t: " << t_max << ", max k: " << k_max << std::endl;
#endif

#if ALWAYS_RESET
    bitset initial_disabled { disabled };
    bitset initial_solved { game.getSolved() };
#endif
    
    while (!pq.empty()) {
        // Step 1: Get values
        auto [k_val, t_val] = pq.top();
        pq.pop();

        // Step 2: Reset the game - we want to know whether this combination can solve the game on its own
        lift_count = 0, lift_attempt = 0;
        uint64_t c;

#if ALWAYS_RESET
        game.reset_to_initial(initial_solved);
        reset_to_initial(initial_disabled);
#endif

#ifndef NDEBUG
        logger << "Currently unsolved: " << game.count_unsolved() << std::endl;
#endif

        // Step 3: Actually do the solving
        if (ODDFIRST) {
            // run odd counters
            run(t_val, k_val, h1, 1);
            c = game.count_unsolved();
#ifndef NDEBUG
            logger << "after odd, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            // if now solved, no need to run odd counters
            if (c != 0)
            {
                // run even counters
                run(t_val, k_val, h0, 0);
                c = game.count_unsolved();
#ifndef NDEBUG
                logger << "after even, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            }
            
        } else {
            // run even counters
            run(t_val, k_val, h0, 0);
            c = game.count_unsolved();
#ifndef NDEBUG
            logger << "after even, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            // if now solved, no need to run odd counters
            if (c != 0)
            {
                // run odd counters
                run(t_val, k_val, h1, 1);
                c = game.count_unsolved();
#ifndef NDEBUG
                logger << "after odd, " << std::setw(9) << lift_count << " lifts, " << std::setw(9) << lift_attempt << " lift attempts, " << c << " unsolved left." << std::endl;
#endif
            }
        }

        // Step 4: Check whether we solved the game
        if (c == 0)
        {
            // We can stop, everything is solved!
            logger << "Solved with k = " << k_val << ", t = " << t_val << std::endl;
            break;
        }
        else if (k_val < k_max or t_val < t_max)
        {
            std::pair<int, int> candidate {k_val + 1, t_val};
            if (k_val + 1 <= k_max and already_tried.find(candidate) == already_tried.end()) 
            {
                pq.push(candidate);
                already_tried.insert(candidate);
            }

            candidate = { k_val, t_val + 1 };
            if (t_val + 1 <= t_max and already_tried.find(candidate) == already_tried.end()) 
            {
                pq.push(candidate);
                already_tried.insert(candidate);
            }
        }
    }
    
}

}
