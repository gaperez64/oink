#!/usr/bin/env python3
"""
Collect NDJSON lines (one JSON object per line) from a file and wrap them
in a JSON array, overwriting the file in-place.

Usage:
    python3 finalize_benchmark.py benchmark_results.json
"""

import json
import sys

def main():
    path = sys.argv[1]
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    with open(path, 'w') as f:
        json.dump(entries, f, indent=2)
    print(f"Wrote {len(entries)} benchmark entries to {path}")

if __name__ == "__main__":
    main()
