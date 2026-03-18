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

// Sentinel values for packed uint8_t level representation.
// TOP_SENTINEL: stored in lane 0 to indicate the "Top" (infinity) progress measure.
// INACTIVE_SENTINEL: stored in unused lanes (lane >= nlanes).
// Valid levels are 0..h-2 (must be < 254 for these sentinels to work).
static constexpr uint8_t TOP_SENTINEL = 0xFF;
static constexpr uint8_t INACTIVE_SENTINEL = 0xFE;

// Bit-parallel popcount for all 8 uint8 lanes simultaneously.
// Uses the standard Hamming-weight algorithm; no intrinsics required.
inline simd_uint8 simd_popcount8(simd_uint8 x) noexcept {
    x = x - ((x >> 1) & simd_uint8(0x55));
    x = (x & simd_uint8(0x33)) + ((x >> 2) & simd_uint8(0x33));
    return (x + (x >> 4)) & simd_uint8(0x0F);
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

    // Flat arrays: pm_bits[node*8 + lane], pm_masks[node*8 + lane], pm_levels[node*8 + lane]
    // This eliminates vector<vector> double-indirection for cache-friendly access.
    std::vector<uint8_t> pm_bits;
    std::vector<uint8_t> pm_masks;
    std::vector<uint8_t> pm_levels;    // uint8_t suffices (levels 0..h-2 < 254)
    std::vector<uint8_t> pm_nlanes;    // number of active lanes per node

    simd_uint8 tmp_bits;
    simd_uint8 tmp_masks;
    simd_uint8 tmp_levels;
    uint8_t tmp_nlanes;

    simd_uint8 best_bits;
    simd_uint8 best_masks;
    simd_uint8 best_levels;
    uint8_t best_nlanes;

    uintqueue Q;
    bitset dirty;

    bool always_reset = false;

    uint64_t *lift_counters;

    // Copy pm[idx] into tmp
    inline void to_tmp(int idx) {
        tmp_bits.copy_from(&pm_bits[idx*8], stdx::element_aligned);
        tmp_masks.copy_from(&pm_masks[idx*8], stdx::element_aligned);
        tmp_levels.copy_from(&pm_levels[idx*8], stdx::element_aligned);
        tmp_nlanes = pm_nlanes[idx];
    }
    // Copy tmp into pm[idx]
    inline void from_tmp(int idx) {
        tmp_bits.copy_to(&pm_bits[idx*8], stdx::element_aligned);
        tmp_masks.copy_to(&pm_masks[idx*8], stdx::element_aligned);
        tmp_levels.copy_to(&pm_levels[idx*8], stdx::element_aligned);
        pm_nlanes[idx] = tmp_nlanes;
    }
    // Copy pm[idx] into best
    inline void to_best(int idx) {
        best_bits.copy_from(&pm_bits[idx*8], stdx::element_aligned);
        best_masks.copy_from(&pm_masks[idx*8], stdx::element_aligned);
        best_levels.copy_from(&pm_levels[idx*8], stdx::element_aligned);
        best_nlanes = pm_nlanes[idx];
    }
    // Copy best into pm[idx]
    inline void from_best(int idx) {
        best_bits.copy_to(&pm_bits[idx*8], stdx::element_aligned);
        best_masks.copy_to(&pm_masks[idx*8], stdx::element_aligned);
        best_levels.copy_to(&pm_levels[idx*8], stdx::element_aligned);
        pm_nlanes[idx] = best_nlanes;
    }
    // Copy tmp into best
    inline void tmp_to_best() {
        best_bits = tmp_bits;
        best_masks = tmp_masks;
        best_levels = tmp_levels;
        best_nlanes = tmp_nlanes;
    }

    // Fill inactive lanes with sentinels
    inline void fill_inactive_tmp() {
        simd_uint8_mask inactive = LANE_INDICES >= simd_uint8(tmp_nlanes);
        stdx::where(inactive, tmp_levels) = simd_uint8(INACTIVE_SENTINEL);
        stdx::where(inactive, tmp_bits) = simd_uint8(0);
        stdx::where(inactive, tmp_masks) = simd_uint8(0);
    }

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render SIMD to given ostream
    void stream_simd(std::ostream &out, simd_uint8& bits, simd_uint8& masks, simd_uint8& levels, uint8_t nlanes);

    // Compare tmp to best
    int compare(int pindex);

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
