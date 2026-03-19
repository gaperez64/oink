window.BENCHMARK_DATA = {
  "lastUpdate": 1773915823129,
  "repoUrl": "https://github.com/gaperez64/oink",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "gaperez64@gmail.com",
            "name": "Guillermo A. Perez",
            "username": "gaperez64"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f536b01b467fc639fe9200f7ea60139381356f03",
          "message": "Claude/add benchmark GitHub action pt igj (#1)\n\n* Add benchmark GitHub Action with regression detection\n\nIntroduces a new CI workflow (.github/workflows/benchmark.yml) that:\n- Builds oink in Release mode on every push/PR to master/main\n- Runs two benchmark suites via test_solvers:\n  * Reproducible random games (fixed seed, 20 games × 5000 nodes) covering\n    fpi, tl, zlk, and pp solver families\n  * A curated subset of the bundled hard-benchmarks (clique, elevator,\n    Friedmann, and real-world CTL/MU-calculus games)\n- Parses the timing output with .github/scripts/parse_benchmark.py and\n  collates results with .github/scripts/finalize_benchmark.py\n- Uses benchmark-action/github-action-benchmark to store historical results\n  in the gh-pages branch and flag regressions (>150% of baseline) with a\n  PR comment and a failing check\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Bump GCC 9 to GCC 11 and add CTest entries for strpm/strpm-simd\n\nGCC 9 predates libstdc++ support for <experimental/simd> (added in\nGCC 11), which causes a compile-time failure for strpm_simd.hpp.\nReplace it with GCC 11 in the CI matrix.\n\nAlso add TestSolverSTRPM and TestSolverSTRPMSIMD to the CTest suite so\nthat both solvers are correctness-verified against the test game set on\nevery CI run.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Focus CTest suite on strpm and strpm-simd only\n\nDisable all other solver tests while developing strpm/strpm-simd to\navoid running the full suite on every ctest invocation.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Restrict benchmark workflow to strpm/strpm-simd; restore all ctest entries\n\nThe regression benchmark now only runs --strpm and --strpm-simd on both\nrandom and hard game sets, keeping CI fast while developing those two\nsolvers.\n\nThe previous commit accidentally dropped all other solver tests from\nCMakeLists.txt; this restores them in full.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Enable CI tests on push/PR; fix benchmark first-run failure\n\n- ci-build.yml: add push/pull_request triggers for master/main so ctest\n  runs automatically (was workflow_dispatch only)\n- benchmark.yml: add a step to initialise the gh-pages branch if absent;\n  the benchmark action errors on first run when the branch does not exist\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Guard strpm_simd behind a feature check for <experimental/simd>\n\nApple libc++ (used by Xcode / macOS CI) does not ship <experimental/simd>;\nonly libstdc++ does.  The macOS jobs in ci-build.yml were silently broken\nbefore because ci-build.yml only ran on workflow_dispatch; now that it\ntriggers on push/PR the failure is visible.\n\n- CMakeLists.txt: use check_include_file_cxx to detect the header and\n  set OINK_HAVE_STRPM_SIMD; only add strpm_simd.cpp and its ctest entry\n  when the header is found.\n- src/solvers.cpp: wrap the strpm_simd include and _add call in\n  #ifdef OINK_HAVE_STRPM_SIMD so the solver is registered only when\n  the implementation was compiled in.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Fix macOS CI: bash syntax, GCC path, stale Xcode entries, checkout@v4\n\nFour pre-existing bugs in the macOS CI that were hidden behind\nworkflow_dispatch but are now exposed with automatic push/PR triggers:\n\n1. Bash syntax error: `if [\"...\"]` → `if [ \"...\" ]` (missing spaces\n   around brackets caused \"command not found\" on every macOS job).\n2. Wrong GCC compiler path: `/usr/bin/gcc` on macOS is a shim to Apple\n   Clang; Homebrew GCC is `gcc-14` / `g++-14` in PATH.\n3. Stale Xcode matrix entries: Xcode 12.0 is unavailable on current\n   macos-latest (macOS 14); no mechanism to switch versions existed\n   anyway. Consolidated to one clang entry (default Xcode).\n4. actions/checkout@v2 uses Node.js 12, which GitHub blocked in 2024;\n   updated to @v4 for both Linux and macOS jobs.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Strengthen SIMD detection: compile-test find_last_set, not just the header\n\ncheck_include_file_cxx only verifies that <experimental/simd> parses.\nstdx::find_last_set() is a GCC libstdc++ extension: Clang on Linux sees\nthe header (via the system GCC libstdc++) but the function may not link\nor compile cleanly.  Replace with check_cxx_source_compiles that tests\nthe exact API the solver relies on, so HAVE_EXPERIMENTAL_SIMD is only\nset when the full required functionality compiles.\n\nAlso switch to cmake_push_check_state/cmake_pop_check_state for clean\nstate management around CMAKE_REQUIRED_FLAGS.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Fix benchmark workflow: git credentials and identity for gh-pages init\n\nTwo issues prevented the gh-pages initialisation step from working:\n\n1. git commit --allow-empty fails without user.name / user.email; add\n   the standard github-actions[bot] identity before any git operations.\n2. git push origin gh-pages needs the token available; pass\n   GITHUB_TOKEN as an env var and use token: in the checkout step\n   (actions/checkout@v4 wires it into the credential helper).\n\nAlso add fetch-depth: 0 so that `git fetch origin gh-pages` can reach\nall remote refs during the existence check.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Fix cmake_push_check_state: include CMakePushCheckState module\n\nThe push/pop check state commands live in CMakePushCheckState.cmake,\nnot CMakeCheckState.cmake.  Without the correct include, cmake\nconfigure failed immediately with \"Unknown CMake command\".\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* Fix gh-pages init: restore HEAD by SHA instead of 'git checkout -'\n\nactions/checkout@v4 leaves the repo in a detached HEAD state, so\n'git checkout -' (return to previous branch) has no branch to go back\nto and fails with 'pathspec - did not match any file(s) known to git'.\n\nSave the HEAD commit SHA before switching to the orphan branch and\nrestore it explicitly afterwards.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n* benchmark CI: enable and run verification tests before benchmarks\n\nEnable OINK_BUILD_TESTS and add a ctest step so the solver correctness\ntests run before the performance benchmarks. This catches regressions in\ncorrectness as well as in performance.\n\nhttps://claude.ai/code/session_01VCprZUoDgs7tBJQ76oK9yA\n\n---------\n\nCo-authored-by: Claude <noreply@anthropic.com>",
          "timestamp": "2026-03-19T11:21:26+01:00",
          "tree_id": "e8546b845238d07811162445a62f7275da020c29",
          "url": "https://github.com/gaperez64/oink/commit/f536b01b467fc639fe9200f7ea60139381356f03"
        },
        "date": 1773915814340,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "strpm [random]",
            "value": 1293,
            "unit": "ms"
          },
          {
            "name": "strpm-simd [random]",
            "value": 1300,
            "unit": "ms"
          },
          {
            "name": "strpm [hard]",
            "value": 94,
            "unit": "ms"
          },
          {
            "name": "strpm-simd [hard]",
            "value": 27860,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "gaperez64@gmail.com",
            "name": "Guillermo A. Perez",
            "username": "gaperez64"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6892cb0c6508d6c320516952fdb59ea887aa479a",
          "message": "Merge branch 'vfluegel:master' into master",
          "timestamp": "2026-03-19T11:21:42+01:00",
          "tree_id": "e8546b845238d07811162445a62f7275da020c29",
          "url": "https://github.com/gaperez64/oink/commit/6892cb0c6508d6c320516952fdb59ea887aa479a"
        },
        "date": 1773915822522,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "strpm [random]",
            "value": 1151,
            "unit": "ms"
          },
          {
            "name": "strpm-simd [random]",
            "value": 1330,
            "unit": "ms"
          },
          {
            "name": "strpm [hard]",
            "value": 49,
            "unit": "ms"
          },
          {
            "name": "strpm-simd [hard]",
            "value": 26855,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}