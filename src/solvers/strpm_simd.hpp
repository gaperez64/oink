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

#ifndef STRPM_SIMD_HPP
#define STRPM_SIMD_HPP

#include "oink/solver.hpp"
#include <experimental/simd>
#include <cstring>

namespace stdx = std::experimental;
using simd_uint8 = stdx::fixed_size_simd<uint8_t, 8>;
using simd_uint8_mask = stdx::fixed_size_simd_mask<uint8_t, 8>;

// Top is represented by levels[0] == -1 (matching the scalar strpm solver).

// Bit-parallel popcount for all 8 uint8 lanes simultaneously.
// Uses the standard Hamming-weight (sideways addition) algorithm:
//   Step 1: Pair adjacent bits: count = b1+b0 for each 2-bit group
//           x - ((x >> 1) & 0x55) works because for 2-bit value ab:
//           ab - 0b = ab (i.e. a+b when both ≤1), handles carry correctly.
//   Step 2: Sum adjacent pairs into 4-bit nibble counts via masking.
//   Step 3: Sum adjacent nibbles; result fits in low nibble, mask off high.
inline simd_uint8 simd_popcount8(simd_uint8 x) noexcept {
    x = x - ((x >> 1) & simd_uint8(0x55));          // 2-bit sums
    x = (x & simd_uint8(0x33)) + ((x >> 2) & simd_uint8(0x33)); // 4-bit sums
    return (x + (x >> 4)) & simd_uint8(0x0F);        // 8-bit sum (0..8)
}

// Precomputed lane index vector [0,1,2,...,7].
static const simd_uint8 LANE_INDICES{[](uint8_t i){ return i; }};

namespace pg {

class STRPM_SIMDSolver : public Solver
{
public:
    STRPM_SIMDSolver(Oink& oink, Game& game);
    virtual ~STRPM_SIMDSolver();

    virtual void run();

protected:
    /**
     * Parameters: U^k_{t, h}
     *      - k: Strahler-number
     *      - t: number of bits
     *      - h: height
     */
    int k, t, h;

    // Flat arrays: pm_bits[node*8 + lane], pm_masks[node*8 + lane]
    // This eliminates vector<vector> double-indirection for cache-friendly access.
    std::vector<uint8_t> pm_bits;
    std::vector<uint8_t> pm_masks;
    // Levels are int (can exceed 255 for large games); flat array stride 8.
    std::vector<int> pm_levels;
    std::vector<uint8_t> pm_nlanes;    // number of active lanes per node

    simd_uint8 tmp_bits;
    simd_uint8 tmp_masks;
    int tmp_levels[8];
    uint8_t tmp_nlanes;

    simd_uint8 best_bits;
    simd_uint8 best_masks;
    int best_levels[8];
    uint8_t best_nlanes;

    uintqueue Q;
    bitset dirty;

    bool always_reset = false;

    uint64_t *lift_counters;

    // Copy pm[idx] into tmp
    inline void to_tmp(int idx) {
        tmp_bits.copy_from(&pm_bits[idx*8], stdx::element_aligned);
        tmp_masks.copy_from(&pm_masks[idx*8], stdx::element_aligned);
        std::memcpy(tmp_levels, &pm_levels[idx*8], 8 * sizeof(int));
        tmp_nlanes = pm_nlanes[idx];
    }
    // Copy tmp into pm[idx]
    inline void from_tmp(int idx) {
        tmp_bits.copy_to(&pm_bits[idx*8], stdx::element_aligned);
        tmp_masks.copy_to(&pm_masks[idx*8], stdx::element_aligned);
        std::memcpy(&pm_levels[idx*8], tmp_levels, 8 * sizeof(int));
        pm_nlanes[idx] = tmp_nlanes;
    }
    // Copy pm[idx] into best
    inline void to_best(int idx) {
        best_bits.copy_from(&pm_bits[idx*8], stdx::element_aligned);
        best_masks.copy_from(&pm_masks[idx*8], stdx::element_aligned);
        std::memcpy(best_levels, &pm_levels[idx*8], 8 * sizeof(int));
        best_nlanes = pm_nlanes[idx];
    }
    // Copy best into pm[idx]
    inline void from_best(int idx) {
        best_bits.copy_to(&pm_bits[idx*8], stdx::element_aligned);
        best_masks.copy_to(&pm_masks[idx*8], stdx::element_aligned);
        std::memcpy(&pm_levels[idx*8], best_levels, 8 * sizeof(int));
        pm_nlanes[idx] = best_nlanes;
    }
    // Copy tmp into best
    inline void tmp_to_best() {
        best_bits = tmp_bits;
        best_masks = tmp_masks;
        std::memcpy(best_levels, tmp_levels, 8 * sizeof(int));
        best_nlanes = tmp_nlanes;
    }

    // Zero out inactive lanes' bits and masks (levels are scalar, not touched)
    inline void fill_inactive_tmp() {
        simd_uint8_mask inactive = LANE_INDICES >= simd_uint8(tmp_nlanes);
        stdx::where(inactive, tmp_bits) = simd_uint8(0);
        stdx::where(inactive, tmp_masks) = simd_uint8(0);
    }

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render SIMD to given ostream
    void stream_simd(std::ostream &out, simd_uint8& bits, simd_uint8& masks, int* levels, uint8_t nlanes);

    // Compare tmp to best (general, works for all k)
    int compare(int pindex);

    // Specialized comparison for k=2: each measure has exactly 1 bitstring,
    // so comparison reduces to scalar operations on single uint8 values.
    // Avoids constructing/comparing full 8-lane SIMD registers.
    // Returns -1 (a < b), 0 (equal), 1 (a > b).
    inline int compare_k2(
        uint8_t a_bits, uint8_t a_masks, int a_level, uint8_t a_nlanes,
        uint8_t b_bits, uint8_t b_masks, int b_level, uint8_t b_nlanes,
        int pindex)
    {
        // Handle Top: represented by nlanes > 0 and level == -1
        if (a_nlanes > 0 and a_level == -1 and b_nlanes > 0 and b_level == -1) return 0;
        if (a_nlanes > 0 and a_level == -1) return 1;
        if (b_nlanes > 0 and b_level == -1) return -1;

        // Handle empty vs non-empty
        if (a_nlanes == 0 and b_nlanes == 0) return 0;
        if (a_nlanes == 0) return (b_bits & 1) ? -1 : 1;
        if (b_nlanes == 0) return (a_bits & 1) ? 1 : -1;

        // Both non-empty, non-Top: compare (level, bitstring) pairs.
        // A string beyond pindex is irrelevant (treated as absent).
        bool a_relevant = (a_level <= pindex);
        bool b_relevant = (b_level <= pindex);
        if (!a_relevant and !b_relevant) return 0;

        // One in range, other out: the one in range has a string, the other doesn't.
        // Leading bit 1 (right subtree) → greater; leading bit 0 → lesser.
        if (a_relevant and !b_relevant) return (a_bits & 1) ? 1 : -1;
        if (!a_relevant and b_relevant) return (b_bits & 1) ? -1 : 1;

        // Both in range: lower level goes first in lexicographic order.
        if (a_level < b_level) return (a_bits & 1) ? 1 : -1;
        if (a_level > b_level) return (b_bits & 1) ? -1 : 1;

        // Same level: full bitstring comparison (masks may differ in length
        // due to partial resets in Case A of prog_tmp).
        // Same algorithm as the SIMD compare(), but on scalar uint8 values.
        uint8_t shorter = a_masks & b_masks;
        uint8_t diff = a_masks ^ b_masks;
        uint8_t first_len_diff = diff & static_cast<uint8_t>(~(diff << 1));
        uint8_t bxor = a_bits ^ b_bits;
        uint8_t combined = shorter + first_len_diff;
        uint8_t rel_xor = bxor & combined;

        // Length-based ordering: shorter string with differing extension
        bool a_less = (a_masks < b_masks and rel_xor == first_len_diff) or
                      (a_masks > b_masks and rel_xor == 0);
        bool b_less = (b_masks < a_masks and rel_xor == first_len_diff) or
                      (b_masks > a_masks and rel_xor == 0);

        // Bit-value ordering within shared prefix
        uint8_t diff_bits = shorter & bxor;
        uint8_t first_bit_diff = diff_bits & (0u - diff_bits); // isolate lowest set bit
        bool a_greater = diff_bits > 0 and (a_bits & first_bit_diff) > 0;
        bool b_greater = diff_bits > 0 and (b_bits & first_bit_diff) > 0;

        if (a_greater or (!a_less and b_less)) return 1;
        if (b_greater or (a_less and !b_less)) return -1;
        return 0;
    }

    void prog_tmp(int pindex, int h);

    // Lift node, triggered by change to target
    bool lift(int node, int target, int &str, int pl);

    inline void todo_push(int node) {
        if (dirty[node]) return;
        Q.push(node);
        dirty[node] = true;
#ifndef NDEBUG
        if (trace >= 2) logger << "push(" << node << ")" << std::endl;
#endif
    }

    inline int todo_pop() {
        int node = Q.pop();
        dirty[node] = false;
#ifndef NDEBUG
        if (trace >= 2) logger << "pop() => " << node << std::endl;
#endif
        return node;
    }

    int lift_count = 0;
    int lift_attempt = 0;

    void run(int t_val, int k_val, int depth, int player);
};

}

#endif
