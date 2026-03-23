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

#include <cstring> // for memset
#include <iomanip>

#include "sspm_simd.hpp"

#define ODDFIRST 1

namespace pg {

SSPM_SIMDSolver::SSPM_SIMDSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

SSPM_SIMDSolver::~SSPM_SIMDSolver()
{
}

/**
 * Set tmp := min { m | m ==_p tmp }
 */
void
SSPM_SIMDSolver::trunc_tmp(int pindex)
{
    if (tmp_levels[0] == -1) return; // already Top
    // compute the lowest pindex >= p
    // [pindex],.,...,.. => [pindex],000
    // if pindex is the bottom, then this simply "buries" the remainder
    int lowest_stp = get_lowest_stp(pindex);
    if (lowest_stp > -1 and lowest_stp < tmp_nlanes - 1) {
        simd_uint16 per_elem = simd_popcount16(tmp_masks);
        // Build inclusive prefix sum: bit_counts[i] = total bits in lanes 0..i
        simd_uint16 bit_counts;
        bit_counts[0] = per_elem[0];
        for (int n = 1; n <= lowest_stp; n++) {
            bit_counts[n] = bit_counts[n-1] + per_elem[n];
        }
        tmp_bits[lowest_stp + 1] = 0;
        // Fill bits so it's exactly l
        int bits_before = lowest_stp >= 0 ? bit_counts[lowest_stp] : 0;
        tmp_masks[lowest_stp + 1] = (1u << (l - bits_before)) - 1;
        tmp_levels[lowest_stp + 1] = pindex + 1;
        tmp_nlanes = lowest_stp + 2;
        fill_inactive_tmp();
    }
}

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
SSPM_SIMDSolver::prog_tmp(int pindex, int h)
{
    // Simple case 1: Top >_p Top
    if (tmp_levels[0] == -1) return; // already Top

    simd_uint16 per_elem = simd_popcount16(tmp_masks);
    // Build inclusive prefix sum: bit_counts[i] = total bits in lanes 0..i
    simd_uint16 bit_counts;
    bit_counts[0] = per_elem[0];
    for (size_t n = 1; n < 16; n++) {
        bit_counts[n] = bit_counts[n-1] + per_elem[n];
    }

    int lowest_stp = get_lowest_stp(pindex);
    // Simple case 2: Some bits below [pindex], ergo [pindex] can go from ..e to ..10*
    if (lowest_stp == -1 or lowest_stp < tmp_nlanes - 1) {
        // We append to an empty string = create a new string
        if (lowest_stp == -1 or tmp_levels[lowest_stp] < pindex)
        {
            lowest_stp ++;
            tmp_bits[lowest_stp] = 0;
            per_elem[lowest_stp] = 0;
            tmp_levels[lowest_stp] = pindex;
        }
        
        // Set the first bit after the current bits to 1
        tmp_bits[lowest_stp] |= (1 << per_elem[lowest_stp]);
        // Fill bits so it's exactly l
        int bits_before = lowest_stp > 0 ? bit_counts[lowest_stp - 1] : 0;
        tmp_masks[lowest_stp] = (1u << (l - bits_before)) - 1;

        tmp_nlanes = lowest_stp + 1; // Add 1 because of 0-indexing
        fill_inactive_tmp();
        return;
    }

    // Case 3: no bits below [pindex], so analyze lowest nonempty level
    // * If lowest contains 0: 3a or 3b
    // * Else if lowest level is root: 3c
    // * Else append 100000000... to next higher level (3d, 3e, 3f)
    //
    // 3a: ,..011*  => ,..100*  (if lowest nonempty is the bottom)
    // 3b: ,..011*, => ,..,000* (if lowest nonempty is not the bottom)
    // 3c: 1111111  => Top      (if root contains only 1s)
    // 3d: ,1111111 => 100*     (if non-root contains only 1s)
    // 3e: ..,111*  => ..100*
    // 3f: ,e,111*  => ,100*
    const auto last_filled = tmp_nlanes-1;
    if (tmp_bits[last_filled] == tmp_masks[last_filled])
    {
        // 3c - f: All 1s
        if (tmp_nlanes == 1)
        {
            if (tmp_levels[0] == 0)
            {
                // 3c: Root with only 1 = go to top
                tmp_bits[0] = 0;
                tmp_masks[0] = 0;
                tmp_levels[0] = -1;
            }
            else
            {
                // 3d: Non-root with only 1s = go up one level and go to 10*
                tmp_bits[0] = 1;
                tmp_levels[0] = tmp_levels[0] - 1;
            }
        }
        else 
        {
            if (tmp_levels[last_filled - 1] + 1 < tmp_levels[last_filled])
            {
                // 3f: there is an empty string to fill = move the string up
                tmp_levels[last_filled] = tmp_levels[last_filled] - 1;
                tmp_bits[last_filled] = 1;
                // Don't have to touch mask - we keep the same number of bits!
            }
            else 
            {
                // 3e: No empty string = append to previous
                tmp_nlanes -= 1;
                // Set the first bit after the current bits to 1
                tmp_bits[last_filled - 1] |= (1 << per_elem[last_filled - 1]);
                // Append the old bits to this mask
                tmp_masks[last_filled - 1] = (1u << (per_elem[last_filled - 1] + per_elem[last_filled])) - 1;
                fill_inactive_tmp();
            }
        }
    }
    else
    {
        // 3 a - b: There is a 0 in there
        // Get the position of the 0
        int last_zero = std::countl_one(static_cast<uint16_t>(tmp_bits[last_filled] | ~tmp_masks[last_filled]));
        if (tmp_levels[last_filled] == h)
        {
            // 3a: Zero in the bottom of the tree = change to 10*
            tmp_bits[last_filled] &= (1u << (15 - last_zero)) - 1;  // Only keep the bits before the 0
            tmp_bits[last_filled] |= (1u << (15 - last_zero));  // Change the 0 to a 1
            // Keep rest as is - Level stays the same, Mask stays
        }
        else 
        {
            // Cut off the existing string at the 0
            tmp_masks[last_filled] = (1u << (15 - last_zero)) - 1;
            int removed_bits = per_elem[last_filled] - std::popcount(static_cast<uint16_t>(tmp_masks[last_filled]));
            tmp_bits[last_filled] &= tmp_masks[last_filled];
            // Add the remaining bits back in the next string as 0s
            tmp_levels[last_filled + 1] = tmp_levels[last_filled] + 1;
            tmp_bits[last_filled + 1] = 0;
            tmp_masks[last_filled + 1] = (1u << removed_bits) - 1;
            tmp_nlanes ++;
        }
    }
}

/**
 * Write pm to ostream.
 */
void
SSPM_SIMDSolver::stream_pm(std::ostream &out, int idx)
{
    // TODO: implement
    /*int base = idx*l;
    if (pm_d[base] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and pm_d[base+j] == i) {
                c++;
                out << pm_b[base+j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }*/
}

/**
 * Write tmp to ostream.
 */
void
SSPM_SIMDSolver::stream_tmp(std::ostream &out, int h)
{
    // TODO: implement
    /*if (tmp_d[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and tmp_d[j] == i) {
                c++;
                out << tmp_b[j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";

        out << " {";

        // compute value
        int i=0;
        for (int d=0; d<h; d++) {
            int val = 0;

            for (; i<l; i++) {
                if (tmp_d[i] != d) {
                    // e found
                    val |= ((1 << (l-i)) - 1);
                    break;
                }

                if (tmp_b[i]) val |= (1 << (l-i));
            }

            logger << " " << val;
        }

        out << " }";
    }*/
}

/**
 * Write best to ostream.
 */
void
SSPM_SIMDSolver::stream_best(std::ostream &out, int h)
{
    // TODO: implement
    /*if (best_d[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and best_d[j] == i) {
                c++;
                out << best_b[j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }*/
}

/**
 * Compare tmp and best
 * res := -1 :: tmp < best
 * res := 0  :: tmp = best
 * res := 1  :: tmp > best
 */
int
SSPM_SIMDSolver::compare(int pindex)
{
    // cases involving Top
    if (tmp_levels[0] == -1 and best_levels[0] == -1) return 0;
    if (tmp_levels[0] == -1) return 1;
    if (best_levels[0] == -1) return -1;

    // Build masks for differences - See strpm_simd for more extensive comments
    simd_uint16 shorter_string = tmp_masks & best_masks;
    simd_uint16 diff = tmp_masks ^ best_masks;
    simd_uint16 first_length_difference = diff & ~(diff << 1);
    simd_uint16 bit_xor = tmp_bits ^ best_bits;
    simd_uint16 combined = shorter_string + first_length_difference;
    simd_uint16 relevant_xor = bit_xor & combined;
    // a_less: tmp's string is shorter (fewer mask bits) and the extra bit in
    //   best is different from tmp (relevant_xor matches the length diff), OR
    //   tmp is longer but all shared bits are identical (the extra 0-extension wins).
    simd_uint16_mask a_less = ((tmp_masks < best_masks) and (relevant_xor == first_length_difference)) or
                             ((tmp_masks > best_masks) and (relevant_xor == 0));
    simd_uint16_mask b_less = ((best_masks < tmp_masks) and (relevant_xor == first_length_difference)) or
                             ((best_masks > tmp_masks) and (relevant_xor == 0));

    simd_uint16 different_bits = (shorter_string & bit_xor);
    simd_uint16 first_bit_difference = different_bits & (simd_uint16(0) - different_bits);
    simd_uint16_mask a_greater = (different_bits > 0) and ((tmp_bits & first_bit_difference) > 0);
    simd_uint16_mask b_greater = (different_bits > 0) and ((best_bits & first_bit_difference) > 0);

    uint16_t max_nlanes = std::max(tmp_nlanes, best_nlanes);
    for (uint16_t i = 0; i < max_nlanes; i++)
    {
        // One of the two has more nonempty strings but we were equal up until here
        if (i >= best_nlanes)
        {
            return (tmp_bits[i] & 1) == 0 ? -1 : 1;
        }
        else if (i >= tmp_nlanes)
        {
            return (best_bits[i] & 1) == 0 ? 1 : -1;
        }
        int tl = tmp_levels[i];
        int bl = best_levels[i];
        // We are past the considered index!
        if (tl > pindex and bl > pindex)
        {
            break;
        }
        // One of the two has a bitstring "earlier"
        else if ((tl <= pindex and bl > pindex) or (tl < bl))
        {
            return (tmp_bits[i] & 1) == 0 ? -1 : 1;
        }
        else if ((tl > pindex and bl <= pindex) or (tl > bl))
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

    return 0;
}

bool
SSPM_SIMDSolver::lift(int v, int target, int &str, int pl)
{
    // check if already Top
    if (pm_levels[16*v] == -1) return false; // already Top

    const int pr = priority(v);
    const int pindex = pl == 0 ? h-(pr+1)/2-1 : h-pr/2-1;

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
                logger << "to target " << label_vertex(target) << ":";
                stream_tmp(logger, h);
                logger << " =>";
            }
#endif
            if (pl == (pr&1)) prog_tmp(pindex, h);
            else trunc_tmp(pindex);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_tmp(logger, h);
                logger << std::endl;
            }
#endif
        to_best(v);
        if (compare(pindex) > 0) {
            from_tmp(v);
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[32;1mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
                stream_tmp(logger, h);
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
            stream_tmp(logger, h);
            logger << " =>";
        }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        else trunc_tmp(pindex);
#ifndef NDEBUG
        if (trace >= 2) {
            stream_tmp(logger, h);
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
            stream_best(logger, h);
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

void
SSPM_SIMDSolver::run(int n_bits, int depth, int player)
{
    l = n_bits;
    h = depth;

    // Initialize progress measures using flat arrays.
    // Every node is set to the smallest leaf in the tree.
    const int nc = nodecount();
    pm_bits.assign(nc * 16, 0);
    pm_levels.assign(nc * 16, 0);
    pm_masks.assign(nc * 16, 0);
    pm_nlanes.assign(nc, static_cast<uint16_t>(1));

    // Initially, every string has all bits in the first level
    uint16_t initial_mask[16] = {};
    initial_mask[0] = (1 << l) - 1;
    for (int n = 0; n < nc; n++)
    {
        std::memcpy(&pm_masks[n * 16], initial_mask, 16);
    }

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
        if (pm_levels[v*16] != -1) {
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

            if (pm_levels[v*16] != -1) {
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
        if (pm_levels[v*16] != -1) Solver::solve(v, 1-player, game.getStrategy(v));
    }

    Solver::flush();

    /*
    delete[] pm_d;
    delete[] tmp_d;
    delete[] best_d;
    delete[] test_d;
    */
}

void
SSPM_SIMDSolver::run()
{
    int max_prio = priority(nodecount()-1);

    // compute ml (max l) and the h for even/odd
    int ml = ceil_log2(nodecount());
    int h0 = (max_prio/2)+1;
    int h1 = (max_prio+1)/2;

    // create datastructures
    Q.resize(nodecount());
    dirty.resize(nodecount());

    logger << "even wants " << ml << "-bounded adaptive " << h0 << "-counters." << std::endl;
    logger << "odd wants " << ml << "-bounded adaptive " << h1 << "-counters." << std::endl;

    // start with 1-bounded adaptive counters
    int i = 1;
    for (; i<=ml; i++) {
        int _l = lift_count, _a = lift_attempt;
        uint64_t _c = game.count_unsolved(), c;

        if (ODDFIRST) {
            // run odd counters
            run(i, h1, 1);
            c = game.count_unsolved();
            logger << "after odd  with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;

            // if now solved, no need to run odd counters
            if (c == 0) break;

            // run even counters
            run(i, h0, 0);
            c = game.count_unsolved();
            logger << "after even with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;
        } else {
            // run even counters
            run(i, h0, 0);
            c = game.count_unsolved();
            logger << "after even with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;

            // if now solved, no need to run odd counters
            if (c == 0) break;

            // run odd counters
            run(i, h1, 1);
            c = game.count_unsolved();
            logger << "after odd  with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;
        }

        if (i != ml) {
            // if i == ml then we are guaranteed to be done
            // otherwise check if done
            if (c == 0) break;
            if (_c != c) i--; // do not increase i if we solved vertices with current i
        } else {
            break; // do not count higher pls
        }
    }

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts, max l " << i << "." << std::endl;
}

}
