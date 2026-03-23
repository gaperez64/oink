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

#ifndef SSPM_SIMD_HPP
#define SSPM_SIMD_HPP

#include "oink/solver.hpp"
#include <experimental/simd>

namespace stdx = std::experimental;
using simd_uint16 = stdx::fixed_size_simd<uint16_t, 16>;
using simd_uint16_mask = stdx::fixed_size_simd_mask<uint16_t, 16>;

namespace pg {

inline simd_uint16 simd_popcount16(simd_uint16 x) noexcept {
    x = x - ((x >> 1) & simd_uint16(0x5555));           // 2-bit sums
    x = (x & simd_uint16(0x3333)) + ((x >> 2) & simd_uint16(0x3333)); // 4-bit sums
    x = (x + (x >> 4)) & simd_uint16(0x0F0F);           // 8-bit sums
    return (x + (x >> 8)) & simd_uint16(0x00FF);        // 16-bit sum (0..16)
}

// Precomputed lane index vector [0,1,2,...,15].
static const simd_uint16 LANE_INDICES_SSPM{[](uint16_t i){ return i; }};

class SSPM_SIMDSolver : public Solver
{
public:
    SSPM_SIMDSolver(Oink& oink, Game& game);
    virtual ~SSPM_SIMDSolver();

    virtual void run();

protected:
    /**
     * So the maximal length of bitstrings l := ceil(log2(n_vertices))
     * And (for even measures) the depth h := # of even priorities
     *
     * So we could reasonably assume 0 <= l <= 32
     *
     * Two formats are possible to encode the array, apart from the bitstring
     * - for each bit the depth: log n * log h 
     */
    int l, h;

    // Flat arrays: pm_bits[node*8 + lane], pm_masks[node*8 + lane]
    // This eliminates vector<vector> double-indirection for cache-friendly access.
    std::vector<uint16_t> pm_bits;
    std::vector<uint16_t> pm_masks;
    // Levels are int (can exceed 255 for large games); flat array stride 8.
    std::vector<int> pm_levels;
    std::vector<uint8_t> pm_nlanes;    // number of active lanes per node

    simd_uint16 tmp_bits;
    simd_uint16 tmp_masks;
    int tmp_levels[16];
    uint8_t tmp_nlanes;

    simd_uint16 best_bits;
    simd_uint16 best_masks;
    int best_levels[16];
    uint8_t best_nlanes;

    uintqueue Q;
    bitset dirty;

    bool bounded = false;

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
        simd_uint16_mask inactive = LANE_INDICES_SSPM >= simd_uint16(tmp_nlanes);
        stdx::where(inactive, tmp_bits) = simd_uint16(0);
        stdx::where(inactive, tmp_masks) = simd_uint16(0);
        for (size_t n = tmp_nlanes; n < 16; n++)
        {
            // Assign max value to all unused levels
            tmp_levels[n] = 0xFFFF;
        }
    }

    // Get the mask for the levels in tmp smaller than p
    inline int get_lowest_stp(int pindex) {
        simd_uint16_mask has_bits = (tmp_masks > 0);
        // smaller than p: check the INT value of the levels against pindex
        alignas(8) uint8_t stp_arr[8] = {};
        for (int i = 0; i < tmp_nlanes; i++) {
            if (tmp_levels[i] <= pindex) stp_arr[i] = 1;
        }
        simd_uint16 stp_v;
        stp_v.copy_from(stp_arr, stdx::element_aligned);
        simd_uint16_mask smaller_than_p = has_bits and (stp_v > simd_uint16(0));
        return stdx::find_last_set(smaller_than_p);
    }

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render tmp to given ostream
    void stream_tmp(std::ostream &out, int h);
    // Render best to given ostream
    void stream_best(std::ostream &out, int h);

    // Compare tmp to best
    int compare(int pindex);

    // Bump tmp given priority p
    void trunc_tmp(int pindex);
    void prog_tmp(int pindex, int h);
    void prog_cap_tmp(int pindex);

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

    void run(int nbits, int depth, int player);
};

}

#endif 
