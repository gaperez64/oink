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

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <queue>

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

// Classify the decisive blocking reason (if any) at the lowest candidate
// coordinate examined by prog_tmp()'s conceptual scalar scan -- i.e. the
// implicit empty coordinate (C1, no SIMD lane) if lowest_is_c1, otherwise
// the deepest explicit lane at or above pindex (lowest_explicit, found via
// find_last_set on smaller_than_p). At most one reason is recorded per
// call; this never alters the masks or the successor result, it is purely
// an observational side channel for scheduling heuristics.
static inline void
record_decisive_block(
    StrpmPressureStats& stats,
    bool lowest_is_c1,
    int lowest_explicit,
    const simd_u16x8_mask& c2,
    const simd_u16x8_mask& c3,
    const simd_u16x8_mask& c4)
{
    if (lowest_is_c1) {
        ++stats.c1_no_string_slot;
        return;
    }
    if (lowest_explicit < 0) return;
    if (c2[lowest_explicit]) {
        ++stats.c2_no_bit_budget;
    } else if (c3[lowest_explicit]) {
        ++stats.c3_zero_ones_boundary;
    } else if (c4[lowest_explicit]) {
        ++stats.c4_ones_boundary;
    }
}

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
STRPM_SIMDSolver::prog_tmp(int pindex, int h)
{
    if (k == 1 and tmp_nlanes == 0)
    {
        if (active_pressure_stats) ++active_pressure_stats->prog_calls;
        // We can immediately handle k = 1, it's just one branch and top
        tmp_levels[0] = TOP_SENTINEL;
        tmp_nlanes = 1;
        tmp_bits = 0;
        tmp_masks = 0;
        fill_inactive_tmp();
        return;
    }

    // Simple case 1: Top >_p Top
    if (tmp_levels[0] == TOP_SENTINEL) return; // already Top
    if (active_pressure_stats) ++active_pressure_stats->prog_calls;

    // NOTE: k=2,t=1 specialization was removed because the assumption that
    // Case A always does a full reset is incorrect: for bit patterns "00" and
    // "01" (bits=0x00/0x01, mask=0x03), countl_one yields 6, giving
    // reset_bits=7 < 8, which triggers a partial reset.

    // --- NLB (Non-Leading-Bit) counting across all 8 bitstring lanes ---
    // Each bitstring is stored as bits/mask pair. The mask has a 1 for each
    // position in the string, and bits holds the actual bit values.
    // NLB count for a single string = popcount(mask) - 1 (subtract the leading bit).
    // For empty strings (mask==0), NLB contribution is 0.
    simd_u16x8_mask has_bits = (tmp_masks > 0);
    simd_u16x8 per_elem = simd_popcount_u16x8(tmp_masks);
    stdx::where(has_bits, per_elem) = per_elem - simd_u16x8(1);
    // Build inclusive prefix sum: nlb_counts[i] = total NLB in lanes 0..i
    // Uses parallel scan (3 dependent steps) on SSE2/NEON, sequential fallback otherwise.
    simd_u16x8 nlb_counts = simd_prefix_sum_inclusive_u16x8(per_elem);
    // nlb_before[i] = total NLB in lanes 0..i-1 (exclusive prefix sum)
    simd_u16x8 nlb_before = nlb_counts - per_elem;

    // --- Per-lane predicates (unified uint16 — no narrowing needed) ---
    // tmp_levels is simd_u16x8, so predicates use it directly — no reload needed.
    // smaller_than_p: lane i's level is at or before the priority index pindex.
    // all_filled_after: from lane i onward, every remaining level slot is occupied.

    // smaller_than_p[i]: tmp_levels[i] <= pindex
    simd_u16x8_mask stp_mask = (tmp_levels <= simd_u16x8(static_cast<uint16_t>(pindex)));

    // all_filled_after[i]: tmp_levels[i] == (h-1-tmp_nlanes) + i
    //   levels are non-decreasing and cover [levels[0], h-2] with no gaps iff this holds.
    simd_u16x8 expected_afa =
        simd_u16x8(static_cast<uint16_t>(h - 1 - static_cast<int>(tmp_nlanes))) + LANE_INDICES;
    simd_u16x8_mask afa_mask = (tmp_levels == expected_afa);

    // With unified uint16, level predicates are directly usable as masks —
    // no int16→uint8 narrowing conversion needed.
    simd_u16x8_mask smaller_than_p = has_bits and stp_mask;

    // --- Determine which lanes have no valid successor (are "stuck") ---
    // clear_first_bit = 0xFFFE: mask to ignore bit position 0 (the leading bit).
    // pattern_zero_and_ones: what the bits look like if leading bit is 0 and
    //   all NLB positions are 1 (i.e., the string is 01^j).
    simd_u16x8 clear_first_bit(static_cast<uint16_t>(0xFFFE));
    simd_u16x8 pattern_zero_and_ones = tmp_masks & clear_first_bit;
    simd_u16x8 t16(static_cast<uint16_t>(t));

    // Named, mutually exclusive blocking-reason masks -- C2/C3/C4 from the
    // handoff (C1 has no SIMD lane, see lowest_is_c1 below). Their union is
    // exactly the original monolithic no_successor formula: with A = "nlb_before
    // == t16", B = "nlb_counts == t16", C = "01^j boundary", D = "1^j boundary",
    // "A or (B and (C or D))" == "A or (!A and B and (C or D))" (X or (!X and Z)
    // == X or Z), and expanding !A-and-B-and-(C-or-D) into the two disjuncts
    // below (with C taking precedence over D) reproduces that exact union.
    const simd_u16x8_mask c2_mask =
        has_bits and (nlb_before == t16);

    const simd_u16x8_mask c3_mask =
        has_bits
        and !c2_mask
        and (nlb_counts == t16)
        and (tmp_bits == pattern_zero_and_ones)
        and afa_mask;

    const simd_u16x8_mask c4_mask =
        has_bits
        and !c2_mask
        and !c3_mask
        and (nlb_counts == t16)
        and (tmp_bits == tmp_masks);

    // A lane has no successor when:
    //   1) nlb_before == t: all t NLB slots are already used in earlier lanes, OR
    //   2) nlb_counts == t (all t NLB used up to and including this lane) AND either:
    //      a) bits == 01^j pattern AND all levels after this are filled (can't reset
    //         downward because there's nowhere to put new strings), OR
    //      b) bits == mask, i.e. all positions are 1 (string is 1^j, maximum value)
    simd_u16x8_mask no_successor = c2_mask or c3_mask or c4_mask;

    // A lane has a successor if it's before pindex, non-empty, and not stuck.
    simd_u16x8_mask has_successor = smaller_than_p and has_bits and !no_successor;

    // --- Pressure telemetry: classify the lowest candidate coordinate ---
    // C1 has no SIMD lane (empty labels aren't stored): if all k-1 nonempty
    // labels already occur above the empty level pindex, that implicit empty
    // coordinate is blocked by C1. Otherwise classify the deepest explicit
    // lane at or above pindex via the masks above. Heuristic-only; does not
    // affect no_successor/has_successor or the successor computation below.
    if (active_pressure_stats) {
        bool lowest_is_c1 = false;
        if (tmp_nlanes == static_cast<uint8_t>(k - 1)) {
            const int last_level = static_cast<int>(tmp_levels[tmp_nlanes - 1]);
            lowest_is_c1 = last_level < pindex;
        }
        const int lowest_explicit = lowest_is_c1 ? -1 : stdx::find_last_set(smaller_than_p);
        record_decisive_block(*active_pressure_stats, lowest_is_c1, lowest_explicit,
                               c2_mask, c3_mask, c4_mask);
    }

    // --- Find the rightmost (highest-index) lane that can be incremented ---
    // We search from the right so that the resulting measure is the smallest
    // value strictly greater than the current one (lexicographic successor).
    int match = stdx::find_last_set(has_successor);
    if (match == -1)
    {
        // No lane can be incremented. Two sub-cases:
        if (tmp_levels[0] == 0)
        {
            // Already at the maximum for this level structure => overflow to Top
            if (active_pressure_stats) ++active_pressure_stats->overflow_to_top;
            tmp_nlanes = 1;
            tmp_levels[0] = TOP_SENTINEL;
            fill_inactive_tmp();
            return;
        }
        else
        {
            // Shift the first string to an earlier level and reset its bits to "1"
            // (the smallest non-empty bitstring with leading bit 1).
            if (active_pressure_stats) ++active_pressure_stats->carry_without_overflow;
            tmp_bits[0] = 1;
            match = 0;
            tmp_levels[0] = static_cast<uint16_t>(std::min(static_cast<int>(tmp_levels[0]) - 1, pindex)); // min needs int for signedness
        }
    }
    else if (nlb_counts[match] == t)
    {
        // Case A: All t NLB slots are consumed up to this lane.
        // The string ends with a pattern like ...10...01^j (trailing ones after
        // the last zero). We must "carry": erase the trailing 01^j suffix.

        // countl_one on (bits | ~mask) counts consecutive 1s from the MSB of
        // the used portion — these are the trailing ones in the bitstring.
        // +1 includes the zero bit that precedes them (the "carry" position).
        // countl_one on uint16: upper 8 bits of ~mask are 0xFF (masks ≤ 8 bits),
        // so countl_one(uint16) counts 8 extra ones; subtract to get the uint8 result.
        int reset_bits = std::countl_one(static_cast<uint16_t>(tmp_bits[match] | ~tmp_masks[match])) - 8 + 1;
        int current_bits = std::popcount(static_cast<uint16_t>(tmp_masks[match]));
        // Clear the top reset_bits positions in both bits and mask:
        // (1 << (8 - reset_bits)) - 1 produces a mask keeping only the lower bits.
        tmp_bits[match] &= (1u << (8-reset_bits)) - 1;
        tmp_masks[match] &= (1u << (8-reset_bits)) - 1;
        if (reset_bits < 8)
        {
            // Partial reset: some bits remain in this lane. Update NLB count
            // and move to the next lane for the tail fill.
            nlb_counts[match] -= current_bits - std::popcount(static_cast<uint16_t>(tmp_masks[match]));
            match++;
            if (match < k-1)
            {
                tmp_levels[match] = tmp_levels[match-1] + 1;
                tmp_bits[match] = 0;
            }
        }
        else
        {
            // Full reset: entire string in this lane was erased.
            // Keep the lane but move it to the next available level.
            nlb_counts[match] -= (current_bits - 1);
            if ((match + 1 == tmp_nlanes and tmp_levels[match] < h-2) or
                (match + 1 < tmp_nlanes and tmp_levels[match] + 1 < tmp_levels[match + 1]) )
            {
                // There's a gap — place this empty string at the next level.
                tmp_levels[match] = (match + 1 == tmp_nlanes)
                    ? tmp_levels[match] + 1
                    : tmp_levels[match + 1] - 1;
            }
            else tmp_levels[match] = tmp_levels[match] + 1;
        }
    }
    else
    {
        // Case B: There are still unused NLB slots available.
        // We can grow the measure by appending a bit.
        if ((match + 1 < tmp_nlanes and tmp_levels[match] + 1 < tmp_levels[match + 1] and tmp_levels[match] != pindex))
        {
            // There's an empty level between this lane and the next —
            // start a new string there with leading bit 1.
            match++;
            tmp_bits[match] = 1;
            tmp_levels[match] = static_cast<uint16_t>(std::min(pindex, static_cast<int>(tmp_levels[match]) - 1)); // min needs int for signedness
        }
        else
        {
            // Extend the current string by setting the next bit position.
            // first_new = current string length = popcount(mask).
            int first_new = std::popcount(static_cast<uint16_t>(tmp_masks[match]));
            tmp_bits[match] |= (1u << first_new);
        }
    }
    if (match < k-1)
    {
        // --- Fill the tail: set all lanes after 'match' to minimum values ---
        // Compute how many NLB are left to distribute, then set this lane's
        // mask to use exactly (t - bits_before + 1) positions (leading bit + NLBs).
        int bits_before = match > 0 ? nlb_counts[match - 1] : 0;
        tmp_masks[match] = (1u << (t - bits_before + 1)) - 1;

        // Append new single-bit strings at consecutive levels after 'match'.
        // Vectorized via simd<uint16_t, 8>: compute the fill range, blend into tmp_levels.
        const int set_to_level = static_cast<int>(tmp_levels[match]) + 1;
        // max_fill: bounded by remaining lane budget (k-2-match) and height budget (h-2-levels[match]).
        // Derivation: original loop ran while nlanes < k-1 AND level <= h-2.
        const int max_fill = std::max(0, std::min(k - 2 - match,
                                                  h - 2 - static_cast<int>(tmp_levels[match])));
        tmp_nlanes = static_cast<uint8_t>(match + 1 + max_fill);
        if (max_fill > 0) {
            // fill_lev[i] = (set_to_level - (match+1)) + i
            //   => fill_lev[match+1] = set_to_level,  fill_lev[match+2] = set_to_level+1, ...
            // set_to_level >= match+1 because levels are non-decreasing (tmp_levels[match] >= match).
            const simd_u16x8 fill_lev =
                simd_u16x8(static_cast<uint16_t>(set_to_level - (match + 1))) + LANE_INDICES;
            // Blend directly into tmp_levels (simd_u16x8 register — no load/store round-trip).
            const simd_u16x8_mask in_range =
                (LANE_INDICES >= simd_u16x8(static_cast<uint16_t>(match + 1))) and
                (LANE_INDICES <  simd_u16x8(static_cast<uint16_t>(match + 1 + max_fill)));
            stdx::where(in_range, tmp_levels) = fill_lev;
        }
        // SIMD bulk-set: for all lanes in (match, tmp_nlanes), set mask=1, bits=0.
        // Lanes beyond tmp_nlanes are zeroed by fill_inactive_tmp below.
        simd_u16x8_mask after_match_mask   = LANE_INDICES > simd_u16x8(static_cast<uint16_t>(match));
        simd_u16x8_mask before_levels_mask = LANE_INDICES < simd_u16x8(tmp_nlanes);
        stdx::where(before_levels_mask and after_match_mask, tmp_masks) = simd_u16x8(1);
        stdx::where(after_match_mask, tmp_bits) = simd_u16x8(0);

        fill_inactive_tmp();
    }
}

/**
 * Write pm to ostream.
 */
void
STRPM_SIMDSolver::stream_pm(std::ostream &out, int idx)
{
    uint8_t nlanes = pm[idx].nlanes;
    if (nlanes > 0 and pm[idx].levels[0] == TOP_SENTINEL) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int output_level = 0;
        for (int i = 0; i < nlanes; i++) {
            if (i>0) out << ", ";
            while (pm[idx].levels[i] > output_level)
            {
                out << "ε, ";
                output_level++;
            }
            out << std::bitset<8>(pm[idx].bits[i]) << "/" << std::bitset<8>(pm[idx].masks[i]);
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
 * Write SIMD state to ostream.
 */
void
STRPM_SIMDSolver::stream_simd(std::ostream &out, simd_u16x8& bits, simd_u16x8& masks, simd_u16x8& levels, uint8_t nlanes)
{
    if (nlanes > 0 and levels[0] == TOP_SENTINEL) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int output_level = 0;
        for (int i = 0; i < nlanes; i++) {
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
    if (k == 1)
    {
        // It is either empty or Top, so comparing sizes is enough
        if (tmp_nlanes == best_nlanes) return 0;
        else if (tmp_nlanes > 0) return 1;
        else if (best_nlanes > 0) return -1;

    }

    // cases involving Top
    if (tmp_levels[0] == TOP_SENTINEL and best_levels[0] == TOP_SENTINEL) return 0;
    if (tmp_levels[0] == TOP_SENTINEL) return 1;
    if (best_levels[0] == TOP_SENTINEL) return -1;

    // --- k=2 specialization: at most 1 active lane, skip SIMD precompute ---
    if (k == 2)
    {
        // Unequal nlanes: one has a string, the other doesn't
        if (tmp_nlanes == 0 and best_nlanes == 0) return 0;
        if (tmp_nlanes > best_nlanes) return (tmp_bits[0] & 1) == 0 ? -1 : 1;
        if (best_nlanes > tmp_nlanes) return (best_bits[0] & 1) == 0 ? 1 : -1;

        // Both have exactly 1 lane. Compare levels then bitstrings.
        int tl = tmp_levels[0], bl = best_levels[0];
        if (tl > pindex and bl > pindex) return 0;
        if ((tl <= pindex and bl > pindex) or (tl < bl))
            return (tmp_bits[0] & 1) == 0 ? -1 : 1;
        if ((tl > pindex and bl <= pindex) or (tl > bl))
            return (best_bits[0] & 1) == 0 ? 1 : -1;

        // Same level: inline scalar bitstring comparison for lane 0.
        uint16_t tm = tmp_masks[0], bm = best_masks[0];
        uint16_t tb = tmp_bits[0],  bb = best_bits[0];
        if (tm == bm and tb == bb) return 0;

        uint16_t shorter = tm & bm;
        uint16_t d = tm ^ bm;
        uint16_t fld = d & ~(d << 1);
        uint16_t bxor = tb ^ bb;
        uint16_t comb = shorter + fld;
        uint16_t rxor = bxor & comb;

        // Length-based comparison
        bool al = (tm < bm and rxor == fld) or (tm > bm and rxor == 0);
        bool bl2 = (bm < tm and rxor == fld) or (bm > tm and rxor == 0);
        // Bit-based comparison within shared prefix
        uint16_t dbits = shorter & bxor;
        uint16_t fbd = dbits & (uint16_t(0) - dbits);
        bool ag = (dbits > 0) and ((tb & fbd) > 0);
        bool bg = (dbits > 0) and ((bb & fbd) > 0);

        if (ag or (!al and bl2)) return 1;
        if (bg or (al and !bl2)) return -1;
        return 0;
    }

    // --- SIMD precomputation: compare all 8 bitstrings in parallel ---
    // Each bitstring is encoded as (bits, mask) where mask indicates which
    // positions exist. A longer string (more mask bits) is compared by first
    // checking the shared prefix, then the extra bit decides ordering.
    //
    // shorter_string: intersection of masks = positions present in both strings.
    // diff: positions that exist in one string but not the other.
    // first_length_difference: isolate the lowest such position (where lengths diverge).
    //   diff & ~(diff << 1) clears all but the least-significant 1 in each
    //   contiguous run of 1s in diff, giving us the first position where one
    //   string is longer than the other.
    simd_u16x8 shorter_string = tmp_masks & best_masks;
    simd_u16x8 diff = tmp_masks ^ best_masks;
    simd_u16x8 first_length_difference = diff & ~(diff << 1);
    simd_u16x8 bit_xor = tmp_bits ^ best_bits;
    // combined = shared positions + the first extra position.
    // relevant_xor = bit differences within this relevant region.
    simd_u16x8 combined     = shorter_string + first_length_difference;
    simd_u16x8 relevant_xor = bit_xor & combined;
    // a_less: tmp's string is shorter (fewer mask bits) and the extra bit in
    //   best is different from tmp (relevant_xor matches the length diff), OR
    //   tmp is longer but all shared bits are identical (the extra 0-extension wins).
    simd_u16x8_mask a_less = ((tmp_masks < best_masks) and (relevant_xor == first_length_difference)) or
                             ((tmp_masks > best_masks) and (relevant_xor == simd_u16x8(0)));
    simd_u16x8_mask b_less = ((best_masks < tmp_masks) and (relevant_xor == first_length_difference)) or
                             ((best_masks > tmp_masks) and (relevant_xor == simd_u16x8(0)));

    // For strings of equal length (or within the shared prefix), find the
    // first bit position where they differ and check who has the 1.
    // different_bits: positions in the shared region where bits disagree.
    // Isolate the lowest such bit with x & (-x) (two's complement trick).
    simd_u16x8 different_bits = (shorter_string & bit_xor);
    simd_u16x8 first_bit_difference = different_bits & (simd_u16x8(0) - different_bits);
    simd_u16x8_mask a_greater = (different_bits > simd_u16x8(0)) and ((tmp_bits & first_bit_difference) > simd_u16x8(0));
    simd_u16x8_mask b_greater = (different_bits > simd_u16x8(0)) and ((best_bits & first_bit_difference) > simd_u16x8(0));

    // --- Vectorized level + string comparison ---
    // pindex can be negative (when priority exceeds the tree depth).
    // uint16_t wraps -1 to 0xFFFF which would make all levels appear <= pindex.
    // Clamp to 0 and use a flag to force all both_active lanes into both_past.
    const bool pindex_negative = (pindex < 0);
    simd_u16x8 pi(static_cast<uint16_t>(std::max(pindex, 0)));

    simd_u16x8_mask in_tmp  = LANE_INDICES < simd_u16x8(tmp_nlanes);
    simd_u16x8_mask in_best = LANE_INDICES < simd_u16x8(best_nlanes);
    simd_u16x8_mask both_active = in_tmp and in_best;

    simd_u16x8_mask tmp_extra  = in_tmp  and !in_best;
    simd_u16x8_mask best_extra = in_best and !in_tmp;

    // When pindex < 0, every level is "past" the index, so all both_active
    // lanes become both_past and none become tl_earlier/bl_earlier.
    simd_u16x8_mask tl_earlier = pindex_negative
        ? simd_u16x8_mask(false)
        : (both_active and
           (((tmp_levels <= pi) and (best_levels > pi)) or (tmp_levels < best_levels)));
    simd_u16x8_mask bl_earlier = pindex_negative
        ? simd_u16x8_mask(false)
        : (both_active and
           (((best_levels <= pi) and (tmp_levels > pi)) or (best_levels < tmp_levels)));

    simd_u16x8_mask both_past = pindex_negative
        ? both_active
        : (both_active and (tmp_levels > pi) and (best_levels > pi));
    simd_u16x8_mask levels_eq = both_active and !tl_earlier and !bl_earlier and !both_past;

    simd_u16x8_mask tmp_lead1  = (tmp_bits  & simd_u16x8(1)) > simd_u16x8(0);
    simd_u16x8_mask best_lead1 = (best_bits & simd_u16x8(1)) > simd_u16x8(0);

    simd_u16x8_mask tmp_decides  = tmp_extra  or tl_earlier;
    simd_u16x8_mask best_decides = best_extra or bl_earlier;

    simd_u16x8_mask a_wins = (tmp_decides  and tmp_lead1)
                          or (best_decides and !best_lead1)
                          or (levels_eq and (a_greater or (!a_less and b_less)));

    simd_u16x8_mask b_wins = (tmp_decides  and !tmp_lead1)
                          or (best_decides and best_lead1)
                          or (levels_eq and (b_greater or (a_less and !b_less)));

    // find_first_set returns a large value (not -1) when no bits are set in
    // libstdc++, so compare against 8 instead of < 0.
    int cutoff = stdx::find_first_set(both_past);
    simd_u16x8_mask before_cutoff = (cutoff >= 8)
        ? simd_u16x8_mask(true)
        : LANE_INDICES < simd_u16x8(static_cast<uint16_t>(cutoff));
    simd_u16x8_mask deciding = (a_wins or b_wins) and before_cutoff;
    int first = stdx::find_first_set(deciding);
    if (first >= 8) return 0;
    return a_wins[first] ? 1 : -1;
}

bool
STRPM_SIMDSolver::lift(int v, int target, int &str, int pl)
{
    // check if already Top
    if (pm[v].nlanes > 0 and pm[v].levels[0] == TOP_SENTINEL) return false; // already Top

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
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << " =>";
            }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        //else trunc_tmp(pindex);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << std::endl;
            }
#endif
        to_best(v);
        if (compare(pindex) > 0) {
            from_tmp(v);
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[32;1mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
                stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
                logger << " (to " << label_vertex(target) << ")\n";
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // Pre-collect valid (non-disabled) successors to avoid branch
    // mispredictions from the disabled check inside the hot loop.
    // succs is a member vector, pre-reserved to nodecount() in run().
    succs.clear();
    for (auto curedge = outs(v); *curedge != -1; curedge++) {
        int to = *curedge;
        if (!disabled[to]) succs.push_back(to);
    }
    int nsuccs = static_cast<int>(succs.size());

    const bool do_prog = (pl == (pr&1));
    const bool want_max = (owner(v) == pl);

    bool first = true;
    for (int si = 0; si < nsuccs; si++) {
        int to = succs[si];

        to_tmp(to);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "to successor " << label_vertex(to) << " from";
            stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
            logger << " =>";
        }
#endif
        if (do_prog) prog_tmp(pindex, h);
#ifndef NDEBUG
        if (trace >= 2) {
            stream_simd(logger, tmp_bits, tmp_masks, tmp_levels, tmp_nlanes);
            logger << std::endl;
        }
#endif
        if (first) {
            tmp_to_best();
            str = to;
            // Early exit: if first successor already at Top and we want max, done
            if (want_max and best_nlanes > 0 and best_levels[0] == TOP_SENTINEL) break;
        } else if (want_max) {
            // Early exit: Top is unsurpassable for max
            if (tmp_nlanes > 0 and tmp_levels[0] == TOP_SENTINEL) {
                tmp_to_best();
                str = to;
                break;
            }
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
            stream_simd(logger, best_bits, best_masks, best_levels, best_nlanes);
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

double
log_binom(int n, int r)
{
    if (r < 0 || r > n) return std::numeric_limits<double>::infinity();
    return std::lgamma(static_cast<double>(n) + 1.0)
         - std::lgamma(static_cast<double>(r) + 1.0)
         - std::lgamma(static_cast<double>(n - r) + 1.0);
}

// Robust log-space size estimate for the Strahler-universal tree:
//   2^(k+t) * binom(t+k-2,k-2) * binom(h-1,k-1)
// Used (neutrally, without pressure weighting) as the base candidate score;
// lower is "smaller tree, cheaper to try". Deliberately does not reuse the
// old constexpr binom()/ApproxSizeCompare -- both had overflow/small-parameter
// issues on this exact formula.
double
log_tree_size_estimate(int k, int t, int h)
{
    if (k <= 1) return 0.0;

    return static_cast<double>(k + t) * std::log(2.0)
         + log_binom(t + k - 2, k - 2)
         + log_binom(h - 1, k - 1);
}

// Is (k,t) inside this schedule's representable/height-bounded grid?
bool
valid(const StrpmOrientationSchedule& schedule, StrpmParams params)
{
    return params.k >= 1 and params.k <= schedule.k_max
       and params.t >= 1 and params.t <= schedule.t_max;
}

// Lazy decrease-key insertion: skip already-tried pairs and non-improving
// scores, otherwise record the new best score for (k,t) and push a fresh
// heap entry (stale entries for the same cell are discarded lazily on pop).
void
enqueue(StrpmOrientationSchedule& schedule, StrpmParams params, double score)
{
    if (!valid(schedule, params)) return;
    if (schedule.tried[params.k][params.t]) return;
    if (score >= schedule.best_pending_score[params.k][params.t]) return;

    schedule.best_pending_score[params.k][params.t] = score;
    schedule.pending.push(StrpmCandidate{params, score, schedule.next_serial++});
}

// Pop the best still-live candidate, discarding stale heap entries whose
// score no longer matches the grid's recorded best for that cell (i.e. a
// later, better enqueue() superseded them) or that were already attempted.
std::optional<StrpmCandidate>
pop_next(StrpmOrientationSchedule& schedule)
{
    while (!schedule.pending.empty()) {
        StrpmCandidate candidate = schedule.pending.top();
        schedule.pending.pop();

        const auto& [k, t] = candidate.params;
        if (schedule.tried[k][t]) continue;
        if (candidate.score != schedule.best_pending_score[k][t]) continue;

        schedule.tried[k][t] = true;
        return candidate;
    }
    return std::nullopt;
}

StrpmDirectionProbabilities
direction_probabilities(const StrpmPressureStats& pressure)
{
    const double k_evidence = static_cast<double>(pressure.k_pressure());
    const double t_evidence = static_cast<double>(pressure.t_pressure());
    const double denominator = k_evidence + t_evidence + 2.0;

    return StrpmDirectionProbabilities{
        (k_evidence + 1.0) / denominator,
        (t_evidence + 1.0) / denominator
    };
}

// After an attempt at (k,t), enqueue both legal neighbors (k+1,t) and
// (k,t+1). Score = estimated log tree size - log(probability that this
// direction addresses the observed blocking pressure); lower is better.
// With no decisive blocks both probabilities are exactly 1/2, so the score
// falls back to the plain tree-size estimate. Both neighbors are always
// offered when in bounds -- pressure only reorders the heap, it never
// prunes a legal direction.
void
expand_after_attempt(StrpmOrientationSchedule& schedule, const StrpmAttemptResult& result)
{
    const int k = result.params.k;
    const int t = result.params.t;

    const StrpmDirectionProbabilities probs = direction_probabilities(result.pressure);

    const double score_k = log_tree_size_estimate(k + 1, t, schedule.h) - std::log(probs.p_increase_k);
    const double score_t = log_tree_size_estimate(k, t + 1, schedule.h) - std::log(probs.p_increase_t);

    if (k + 1 <= schedule.k_max)
        enqueue(schedule, {k + 1, t}, score_k);

    if (t + 1 <= schedule.t_max)
        enqueue(schedule, {k, t + 1}, score_t);
}

StrpmAttemptResult
STRPM_SIMDSolver::run_attempt(int t_val, int k_val, int depth, int player)
{
    StrpmAttemptResult result;
    result.player = player;
    result.params = {k_val, t_val};
    result.h = depth + 1;
    result.unsolved_before = static_cast<uint64_t>(game.count_unsolved());

    // Marcin's word: think of h as the number of priorities of the
    // opponent... PLUS ONE!
    t = t_val;
    h = depth + 1;  // FIXME: This is Guillermo's hack, the +1
    k = k_val;  // Maybe possible: std::min(t + 2, h);

    lift_count = 0;
    lift_attempt = 0;

#ifndef NDEBUG
    logger << "Strahler-tree parameters for player " << player << ": k = " << k << ", t = " << t << ", h = " << h << std::endl;
#endif

    // Initialize progress measures using AoS layout.
    // Every node is set to the smallest leaf in the tree.
    const int nc = nodecount();

    // Build the initial NodePM pattern once, then fill the entire vector with it.
    NodePM init{};
    init.nlanes  = static_cast<uint8_t>(k - 1);
    init.masks[0] = static_cast<uint16_t>((1 << (t+1)) - 1);
    // init.levels[0] = 0 already (zero-initialised by NodePM{})
    for (int i = 1; i < k-1; i++)
    {
        init.levels[i] = static_cast<uint16_t>(i);
        init.masks[i]  = 1;
    }
    pm.assign(nc, init);

#ifndef NDEBUG
    if (trace >= 1)
    {
        logger << "Initial PM: " << std::endl;
        stream_pm(logger, 0);
        logger << std::endl;
    }
#endif

    // Pressure telemetry covers only the initial-lifting fixed point below;
    // the strategy-extraction re-lift further down must not pollute it (see
    // the active_pressure_stats = nullptr reset before that pass).
    active_pressure_stats = &result.pressure;

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

    // This re-lifts every vertex to extract strategies, not to explore the
    // fixed point -- it must not be counted as scheduling pressure.
    active_pressure_stats = nullptr;

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm[v].nlanes == 0 or pm[v].levels[0] != TOP_SENTINEL) {
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

            if (pm[v].nlanes == 0 or pm[v].levels[0] != TOP_SENTINEL) {
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
        if (pm[v].nlanes == 0 or pm[v].levels[0] != TOP_SENTINEL) Solver::solve(v, 1-player, game.getStrategy(v));
    }

    Solver::flush();

    result.unsolved_after = static_cast<uint64_t>(game.count_unsolved());
    result.lifts = static_cast<uint64_t>(lift_count);
    result.lift_attempts = static_cast<uint64_t>(lift_attempt);
    return result;
}

void
STRPM_SIMDSolver::run()
{
    int max_prio = priority(nodecount()-1);

    // compute the h for even/odd
    int h0 = (max_prio/2)+1;
    int h1 = (max_prio+1)/2;

    assert(std::max(h0, h1) < 65535); // levels stored as uint16_t; 0xFFFF reserved for Top

    // create datastructures
    Q.resize(nodecount());
    dirty.resize(nodecount());
    succs.reserve(nodecount());

    // Two independent per-orientation schedules replace the old single
    // global (k,t) queue: each has its own height, candidate heap,
    // attempted/pending grids, and pressure history (via last_result), and
    // a run for player p only ever touches schedules[p].
    std::array<StrpmOrientationSchedule, 2> schedules;
    schedules[0].player = 0;
    schedules[0].h = h0 + 1; // the exact h run_attempt(..., h0, ...) uses (depth+1)
    schedules[1].player = 1;
    schedules[1].h = h1 + 1;

    for (auto& schedule : schedules) {
        // Representable region of the SIMD encoding. Even with the uint16
        // unification, bitstrings are hard-capped at 8 bits (the mask/countl_one
        // logic uses a literal 8), so t <= STRPM_SIMD_T_MAX; there are 8 SIMD
        // lanes and nlanes == k-1, so k <= STRPM_SIMD_K_MAX. Beyond this the
        // encoding would silently overflow, so we cap the (k,t) search here
        // (per orientation, bounded by that orientation's own height) and
        // hand any still-unsolved remainder to a complete solver (tangle
        // learning) after the loop.
        schedule.t_max = std::min(floor_log2(nodecount()), STRPM_SIMD_T_MAX);
        schedule.k_max = std::min({schedule.t_max + 2, schedule.h, STRPM_SIMD_K_MAX});
        // Same starting pair the old global schedule used.
        enqueue(schedule, {1, 1}, log_tree_size_estimate(1, 1, schedule.h));
    }

#ifndef NDEBUG
    logger << "player 0: max t: " << schedules[0].t_max << ", max k: " << schedules[0].k_max << std::endl;
    logger << "player 1: max t: " << schedules[1].t_max << ", max k: " << schedules[1].k_max << std::endl;
#endif

#if ALWAYS_RESET
    bitset initial_disabled { disabled };
    bitset initial_solved { game.getSolved() };
#endif

    // Fair round-robin dispatcher: each round performs at most one attempt
    // per orientation, so neither orientation can starve the other by
    // generating an endless chain of favored candidates. The schedules
    // remain valid as the other orientation solves vertices, since the
    // residual game only shrinks.
    const std::array<int, 2> order =
        ODDFIRST ? std::array<int,2>{1,0}
                 : std::array<int,2>{0,1};

    while (game.count_unsolved() != 0) {
        bool ran_any_attempt = false;

        for (int player : order) {
            auto& schedule = schedules[player];
            const auto candidate = pop_next(schedule);
            if (!candidate) continue;

            ran_any_attempt = true;

#if ALWAYS_RESET
            game.reset_to_initial(initial_solved);
            reset_to_initial(initial_disabled);
#endif

            const auto result = run_attempt(
                candidate->params.t,
                candidate->params.k,
                player == 0 ? h0 : h1,
                player);

            schedule.last_result = result;
            expand_after_attempt(schedule, result);

#ifndef NDEBUG
            // One compact summary line per orientation attempt; never
            // per-call prog_tmp() pressure detail (that would be too hot).
            if (trace >= 1) {
                const int k = candidate->params.k;
                const int t = candidate->params.t;
                logger << "strpm-simd p=" << player << " k=" << k << " t=" << t
                       << " h=" << result.h << std::endl;
                logger << "  solved=" << result.solved_vertices()
                       << " remaining=" << result.unsolved_after << std::endl;
                logger << "  lifts=" << result.lifts << " attempts=" << result.lift_attempts << std::endl;
                logger << "  pressure={c1:" << result.pressure.c1_no_string_slot
                       << ",c2:" << result.pressure.c2_no_bit_budget
                       << ",c3:" << result.pressure.c3_zero_ones_boundary
                       << ",c4:" << result.pressure.c4_ones_boundary
                       << ",top:" << result.pressure.overflow_to_top << "}" << std::endl;

                const StrpmDirectionProbabilities probs = direction_probabilities(result.pressure);
                if (k + 1 <= schedule.k_max) {
                    logger << "  enqueue k+1 score="
                           << (log_tree_size_estimate(k + 1, t, schedule.h) - std::log(probs.p_increase_k))
                           << std::endl;
                }
                if (t + 1 <= schedule.t_max) {
                    logger << "  enqueue t+1 score="
                           << (log_tree_size_estimate(k, t + 1, schedule.h) - std::log(probs.p_increase_t))
                           << std::endl;
                }
            }
#endif

            if (game.count_unsolved() == 0) {
                logger << "Solved with k = " << candidate->params.k << ", t = " << candidate->params.t << std::endl;
                break;
            }
        }

        if (!ran_any_attempt)
            break;
    }

    // Fallback: each schedule's (k,t) search is capped at the representable
    // region (k <= STRPM_SIMD_K_MAX, t <= STRPM_SIMD_T_MAX) and its own
    // height-bounded k_max/t_max. If both schedules drained without solving
    // the whole game, hand the remaining unsolved subgame to a complete
    // solver. All vertices strpm-simd already solved are in <disabled>, so
    // tangle learning only finishes what is left.
    if (game.count_unsolved() != 0)
    {
        for (int player = 0; player < 2; player++) {
            const auto& schedule = schedules[player];
            logger << "strpm-simd exhausted player " << player << " schedule up to k = "
                   << schedule.last_result.params.k << ", t = " << schedule.last_result.params.t
                   << " (bounds k <= " << schedule.k_max << ", t <= " << schedule.t_max << ")" << std::endl;
        }
        logger << game.count_unsolved()
               << " vertices unsolved; falling back to tangle learning." << std::endl;
        solveRemainderWith("tl");
    }
}

}
