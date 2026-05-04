#!/bin/bash
# bfqutils test suite orchestrator.
# Run from the test/ directory, or via `make test` from the repo root.

cd "$(dirname "$0")" || exit 1
source lib.sh

total_pass=0
total_fail=0

run_suite() {
    local script=$1
    PASS=0 FAIL=0
    source "$script"
    total_pass=$((total_pass + PASS))
    total_fail=$((total_fail + FAIL))
}

run_suite tests/test_merge.sh
run_suite tests/test_trimse.sh
run_suite tests/test_trimpe.sh
run_suite tests/test_stats.sh

echo ""
echo "=== Total: $total_pass/$((total_pass + total_fail)) passed ==="
[[ $total_fail -eq 0 ]]
