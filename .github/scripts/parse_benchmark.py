#!/usr/bin/env python3
"""
Parse test_solvers output and emit benchmark JSON objects for
github-action-benchmark (customSmallerIsBetter format).

Usage:
    ./build/test_solvers ... | python3 parse_benchmark.py <label>

Each invocation appends one JSON object per solver to stdout.
The label is appended to the benchmark name, e.g. "fpi [random]".
"""

import re
import sys
import json

ANSI_ESCAPE = re.compile(r'\x1b\[[0-9;]*m')

def strip_ansi(s: str) -> str:
    return ANSI_ESCAPE.sub('', s)

def parse_times_line(line: str):
    """
    Parse a line like:
        times:   fpi (1234 ms) tl (567 ms) zlk (89 ms)
    Returns a list of (solver_name, ms_value) tuples.
    """
    # Remove 'times:' prefix
    line = re.sub(r'^times:\s*', '', line)
    # Find all  <name> (<value> ms)  patterns
    return re.findall(r'(\S+)\s+\((\d+(?:\.\d+)?)\s+ms\)', line)

def main():
    label = sys.argv[1] if len(sys.argv) > 1 else "benchmark"
    results = []

    for raw_line in sys.stdin:
        clean = strip_ansi(raw_line)
        if clean.startswith('times:'):
            for solver, ms in parse_times_line(clean):
                results.append({
                    "name": f"{solver} [{label}]",
                    "value": float(ms),
                    "unit": "ms"
                })
            break  # only one times: line per run

    for entry in results:
        print(json.dumps(entry))

if __name__ == "__main__":
    main()
