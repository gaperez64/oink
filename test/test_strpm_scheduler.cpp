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

// Pure-function tests for the strpm-simd player-specific, failure-directed
// scheduler (log_binom/log_tree_size_estimate, direction_probabilities,
// enqueue/pop_next/expand_after_attempt, and per-orientation schedule
// independence). No game or solver is constructed -- this exercises only
// the scheduling helpers declared in solvers/strpm_simd.hpp. Mirrors the
// hand-rolled harness style of test_solvers.cpp: no external test
// framework, plain checks with a running failure count and a nonzero exit
// code on failure.

#include "solvers/strpm_simd.hpp"

#include <cmath>
#include <iostream>

using namespace pg;

namespace {

int failures = 0;

void
check(bool cond, const char* what)
{
    if (!cond) {
        std::cerr << "FAIL: " << what << std::endl;
        failures++;
    }
}

StrpmPressureStats
make_pressure(uint64_t c1, uint64_t c2, uint64_t c3, uint64_t c4)
{
    StrpmPressureStats p;
    p.c1_no_string_slot = c1;
    p.c2_no_bit_budget = c2;
    p.c3_zero_ones_boundary = c3;
    p.c4_ones_boundary = c4;
    return p;
}

StrpmAttemptResult
make_result(int k, int t, StrpmPressureStats pressure)
{
    StrpmAttemptResult r;
    r.params = {k, t};
    r.pressure = pressure;
    return r;
}

void
test_log_binom_sanity()
{
    check(std::isfinite(log_binom(5, 2)), "log_binom(5,2) is finite");
    check(std::isinf(log_binom(5, 6)), "log_binom(5,6) is out of range -> +inf");
    check(std::isinf(log_binom(5, -1)), "log_binom(5,-1) is out of range -> +inf");
    check(std::abs(log_binom(5, 0) - 0.0) < 1e-9, "log_binom(n,0) == log(1) == 0");
    check(log_binom(10, 5) > log_binom(6, 3), "log_binom grows with n at a fixed ratio");
}

void
test_log_tree_size_estimate_sanity()
{
    check(log_tree_size_estimate(1, 1, 10) == 0.0, "k<=1 collapses to 0");
    check(log_tree_size_estimate(0, 1, 10) == 0.0, "k<=1 collapses to 0 (k=0)");
    check(std::isfinite(log_tree_size_estimate(3, 2, 10)), "finite for in-range params");
    check(log_tree_size_estimate(3, 2, 20) > log_tree_size_estimate(3, 2, 10),
          "estimate grows with height h");
    check(log_tree_size_estimate(5, 2, 20) > log_tree_size_estimate(3, 2, 20),
          "estimate grows with k");
    check(log_tree_size_estimate(3, 5, 20) > log_tree_size_estimate(3, 2, 20),
          "estimate grows with t");
}

void
test_direction_probabilities()
{
    {
        const auto probs = direction_probabilities(make_pressure(0, 0, 0, 0));
        check(std::abs(probs.p_increase_k - 0.5) < 1e-9, "no blockers => p_increase_k == 0.5");
        check(std::abs(probs.p_increase_t - 0.5) < 1e-9, "no blockers => p_increase_t == 0.5");
    }
    {
        const auto probs = direction_probabilities(make_pressure(10, 0, 0, 0));
        check(probs.p_increase_k > probs.p_increase_t, "C1 only favors growing k");
    }
    {
        const auto probs = direction_probabilities(make_pressure(0, 10, 0, 0));
        check(probs.p_increase_t > probs.p_increase_k, "C2 only favors growing t");
    }
    {
        const auto probs = direction_probabilities(make_pressure(0, 0, 10, 0));
        check(probs.p_increase_t > probs.p_increase_k, "C3 only favors growing t");
    }
    {
        const auto probs = direction_probabilities(make_pressure(0, 0, 0, 10));
        check(probs.p_increase_t > probs.p_increase_k, "C4 only favors growing t");
    }
    {
        // Laplace smoothing: neither probability ever reaches exactly 0 or 1,
        // even under extreme one-sided pressure.
        const auto probs = direction_probabilities(make_pressure(1000000, 0, 0, 0));
        check(probs.p_increase_t > 0.0 && probs.p_increase_t < 1.0,
              "p_increase_t stays in (0,1) under heavy one-sided C1 pressure");
        check(probs.p_increase_k > 0.0 && probs.p_increase_k < 1.0,
              "p_increase_k stays in (0,1) under heavy one-sided C1 pressure");
    }
}

void
test_candidate_insertion_completeness()
{
    // Both neighbors pending when both (k+1,t) and (k,t+1) are in bounds.
    {
        StrpmOrientationSchedule sched;
        sched.h = 20;
        sched.k_max = 9;
        sched.t_max = 7;
        expand_after_attempt(sched, make_result(2, 2, make_pressure(0, 0, 0, 0)));
        check(sched.pending.size() == 2, "both neighbors enqueued when both in bounds");
    }
    // Only the t-neighbor when k is already at its cap.
    {
        StrpmOrientationSchedule sched;
        sched.h = 20;
        sched.k_max = 2;
        sched.t_max = 7;
        expand_after_attempt(sched, make_result(2, 2, make_pressure(0, 0, 0, 0)));
        const auto c = pop_next(sched);
        check(c.has_value() && c->params.k == 2 && c->params.t == 3,
              "only the t-neighbor is enqueued when k is capped");
        check(!pop_next(sched).has_value(), "no second candidate when k was capped");
    }
    // Only the k-neighbor when t is already at its cap.
    {
        StrpmOrientationSchedule sched;
        sched.h = 20;
        sched.k_max = 9;
        sched.t_max = 2;
        expand_after_attempt(sched, make_result(2, 2, make_pressure(0, 0, 0, 0)));
        const auto c = pop_next(sched);
        check(c.has_value() && c->params.k == 3 && c->params.t == 2,
              "only the k-neighbor is enqueued when t is capped");
        check(!pop_next(sched).has_value(), "no second candidate when t was capped");
    }
    // At k=1, prog_tmp()'s k==1 fast path ignores t entirely (every measure
    // is either the initial empty value or Top), so (1,t) and (1,t+1) are
    // provably the same attempt. Only the k-neighbor should be offered.
    {
        StrpmOrientationSchedule sched;
        sched.h = 20;
        sched.k_max = 9;
        sched.t_max = 7;
        expand_after_attempt(sched, make_result(1, 1, make_pressure(0, 0, 0, 0)));
        const auto c = pop_next(sched);
        check(c.has_value() && c->params.k == 2 && c->params.t == 1,
              "from k=1, only the k-neighbor is enqueued (t-growth is a no-op at k=1)");
        check(!pop_next(sched).has_value(),
              "no (1,t+1) candidate is ever offered from a k=1 attempt");
    }
    // No pair is ever attempted twice.
    {
        StrpmOrientationSchedule sched;
        sched.h = 20;
        sched.k_max = 9;
        sched.t_max = 7;
        enqueue(sched, {1, 1}, 0.0);
        const auto c1 = pop_next(sched);
        check(c1.has_value() && c1->params.k == 1 && c1->params.t == 1, "first pop returns (1,1)");
        enqueue(sched, {1, 1}, -100.0); // re-offer an already-tried pair with a far better score
        check(!pop_next(sched).has_value(), "an already-tried pair is never re-attempted");
    }
    // A better later score creates a new heap entry; the stale old entry is
    // discarded lazily (on pop), not eagerly.
    {
        StrpmOrientationSchedule sched;
        sched.h = 20;
        sched.k_max = 9;
        sched.t_max = 7;
        enqueue(sched, {2, 2}, 10.0);
        enqueue(sched, {2, 2}, 5.0);  // better -- supersedes 10.0
        enqueue(sched, {2, 2}, 8.0);  // worse than the current best (5.0) -- ignored
        check(sched.pending.size() == 2, "two heap entries pushed (10.0 stale once superseded, 5.0 live)");

        const auto c = pop_next(sched);
        check(c.has_value() && std::abs(c->score - 5.0) < 1e-9,
              "pop_next returns the better (later) score, not the first-pushed one");
        check(sched.pending.size() == 1, "the stale 10.0 entry is still in the heap after popping 5.0");
        check(!pop_next(sched).has_value(),
              "the stale 10.0 entry is discarded lazily on its own pop, not returned");
        check(sched.pending.empty(), "heap is empty once the stale entry has been discarded");
    }
}

void
test_schedule_independence()
{
    StrpmOrientationSchedule even_sched;
    even_sched.player = 0;
    even_sched.h = 20;
    even_sched.k_max = 9;
    even_sched.t_max = 7;

    StrpmOrientationSchedule odd_sched;
    odd_sched.player = 1;
    odd_sched.h = 25;
    odd_sched.k_max = 9;
    odd_sched.t_max = 7;

    enqueue(odd_sched, {1, 1}, 0.0);
    const auto odd_size_before = odd_sched.pending.size();
    const auto odd_serial_before = odd_sched.next_serial;
    const bool odd_tried_2_2_before = odd_sched.tried[2][2];

    expand_after_attempt(even_sched, make_result(2, 2, make_pressure(5, 0, 0, 0)));

    check(odd_sched.pending.size() == odd_size_before,
          "expanding the even schedule does not change the odd schedule's heap size");
    check(odd_sched.next_serial == odd_serial_before,
          "expanding the even schedule does not change the odd schedule's serial counter");
    check(odd_sched.tried[2][2] == odd_tried_2_2_before,
          "expanding the even schedule does not change the odd schedule's tried grid");
    check(even_sched.pending.size() == 2, "the even schedule got its own two neighbors");

    // Different pressure profiles can steer two schedules toward different
    // next pairs from the same starting (k,t).
    StrpmOrientationSchedule k_favored;
    k_favored.h = 20;
    k_favored.k_max = 9;
    k_favored.t_max = 7;
    expand_after_attempt(k_favored, make_result(2, 2, make_pressure(1000, 0, 0, 0)));

    StrpmOrientationSchedule t_favored;
    t_favored.h = 20;
    t_favored.k_max = 9;
    t_favored.t_max = 7;
    expand_after_attempt(t_favored, make_result(2, 2, make_pressure(0, 1000, 0, 0)));

    const auto next_k = pop_next(k_favored);
    const auto next_t = pop_next(t_favored);
    check(next_k.has_value() && next_k->params.k == 3 && next_k->params.t == 2,
          "heavy C1 pressure steers the schedule to try (k+1,t) first");
    check(next_t.has_value() && next_t->params.k == 2 && next_t->params.t == 3,
          "heavy C2 pressure steers the schedule to try (k,t+1) first");
}

} // namespace

int
main()
{
    test_log_binom_sanity();
    test_log_tree_size_estimate_sanity();
    test_direction_probabilities();
    test_candidate_insertion_completeness();
    test_schedule_independence();

    if (failures == 0) {
        std::cout << "test_strpm_scheduler: all checks passed" << std::endl;
        return 0;
    }
    std::cerr << "test_strpm_scheduler: " << failures << " check(s) failed" << std::endl;
    return 1;
}
