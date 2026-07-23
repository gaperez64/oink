#!/bin/bash
# strpm-simd scheduling-diagnostics tests (handoff section 11.3): four tiny
# deterministic games, each expected to make at least one attempt report a
# specific blocking-pressure category (or none). For each fixture:
#   - the solver must run to completion and Oink's Verifier must accept the
#     resulting strategy (the project's sole correctness oracle, matching
#     test_solvers.cpp's convention);
#   - the aggregate counter for the fixture's named category, summed across
#     every attempt's trace>=1 "pressure={...}" line, must be nonzero (or,
#     for the "no pressure" fixture, every category must be zero).
# Exact counts are deliberately not asserted -- only nonzero-ness -- since
# scheduling order (and therefore exact per-attempt counts) is not part of
# the correctness contract.
set -u

OINK_BIN="$1"
FIXTURE_DIR="$2"

failures=0

sum_counter() {
    # $1 = trace output, $2 = counter name (c1|c2|c3|c4)
    echo "$1" | grep -oE "$2:[0-9]+" | cut -d: -f2 | awk '{s+=$1} END {print s+0}'
}

check_fixture() {
    local file="$1" label="$2"
    local out
    out=$("$OINK_BIN" --strpm-simd -t --no -v "$file" 2>&1)
    local rc=$?

    if [ $rc -ne 0 ] || ! echo "$out" | grep -qi "solution verified"; then
        echo "FAIL [$label]: solver did not run to completion / verifier did not accept the solution (rc=$rc)"
        echo "$out" | tail -10
        failures=$((failures+1))
        return
    fi

    local c1 c2 c3 c4
    c1=$(sum_counter "$out" c1)
    c2=$(sum_counter "$out" c2)
    c3=$(sum_counter "$out" c3)
    c4=$(sum_counter "$out" c4)
    echo "  [$label] aggregate pressure: c1=$c1 c2=$c2 c3=$c3 c4=$c4"

    case "$label" in
        c1)
            [ "$c1" -gt 0 ] || { echo "FAIL [$label]: expected nonzero c1_no_string_slot"; failures=$((failures+1)); }
            ;;
        c2)
            [ "$c2" -gt 0 ] || { echo "FAIL [$label]: expected nonzero c2_no_bit_budget"; failures=$((failures+1)); }
            ;;
        c3_or_c4)
            [ "$((c3 + c4))" -gt 0 ] || { echo "FAIL [$label]: expected nonzero c3_zero_ones_boundary or c4_ones_boundary"; failures=$((failures+1)); }
            ;;
        none)
            [ "$((c1 + c2 + c3 + c4))" -eq 0 ] || { echo "FAIL [$label]: expected no decisive blocking pressure at all"; failures=$((failures+1)); }
            ;;
    esac
}

check_fixture "$FIXTURE_DIR/c1_pressure.pg"     c1
check_fixture "$FIXTURE_DIR/c2_pressure.pg"     c2
check_fixture "$FIXTURE_DIR/c3_c4_pressure.pg"  c3_or_c4
check_fixture "$FIXTURE_DIR/no_pressure.pg"     none

if [ "$failures" -eq 0 ]; then
    echo "test_strpm_pressure: all fixtures matched their expected pressure profile"
    exit 0
fi
echo "test_strpm_pressure: $failures fixture(s) failed"
exit 1
