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

namespace stdx = std::experimental;
using simd_uint8 = stdx::fixed_size_simd<uint8_t, 8>;
using simd_uint8_mask = stdx::fixed_size_simd_mask<uint8_t, 8>;

// Bit-parallel popcount for all 8 uint8 lanes simultaneously.
// Uses the standard Hamming-weight algorithm; no intrinsics required.
inline simd_uint8 simd_popcount8(simd_uint8 x) noexcept {
    x = x - ((x >> 1) & simd_uint8(0x55));
    x = (x & simd_uint8(0x33)) + ((x >> 2) & simd_uint8(0x33));
    return (x + (x >> 4)) & simd_uint8(0x0F);
}

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
    std::vector<std::vector<uint8_t>> pm_bits;
    std::vector<std::vector<uint8_t>> pm_masks;
    std::vector<std::vector<int>> pm_levels;

    simd_uint8 tmp_bits;
    simd_uint8 tmp_masks;
    std::vector<int> tmp_levels;

    simd_uint8 best_bits;
    simd_uint8 best_masks;
    std::vector<int> best_levels;

    uintqueue Q;
    bitset dirty;

    bool always_reset = false;

    uint64_t *lift_counters;

    // Copy pm[idx] into tmp
    void to_tmp(int idx);
    // Copy tmp into pm[idx]
    void from_tmp(int idx);
    // Copy pm[idx] into best
    void to_best(int idx);
    // Copy best into pm[idx]
    void from_best(int idx);
    // Copy tmp into best
    void tmp_to_best();

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render SIMD to given ostream
    void stream_simd(std::ostream &out, simd_uint8& bits, simd_uint8& masks, std::vector<int>& levels);

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
