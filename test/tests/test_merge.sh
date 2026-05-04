#!/bin/bash
# Tests for bfqutils merge
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source lib.sh
BFQ=$BFQUTILS

echo "=== merge ==="

# Regression: default merge output unchanged
run_diff "default merge" \
    expected/merge_default.fq \
    $BFQ merge -q fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Tighter overlap requirement produces fewer merged reads; omit -q to see "Created" summary
run_grep "overlap -o 25 succeeds" "Created" \
    $BFQ merge -o 25 fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# -g flag reaches the polyG code (was dead code before v1.2.2)
run_ok "-g flag accepted" \
    $BFQ merge -q -g 8 fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Comma-separated input: two pairs → output is double the single-pair run
single_n=$($BFQ merge -q fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | wc -l)
double_n=$($BFQ merge -q \
    fixtures/testR1.fq.gz,fixtures/testR1.fq.gz \
    fixtures/testR2.fq.gz,fixtures/testR2.fq.gz 2>/dev/null | wc -l)
if [[ $((single_n * 2)) -eq "$double_n" ]]; then
    _log_pass "comma-separated inputs double output"
else
    _log_fail "comma-separated inputs: expected $((single_n * 2)) lines, got $double_n"
fi

# Gzip output (-z): decompress and diff against plain reference
gz_out=$($BFQ merge -q -z fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | gunzip)
if diff -q expected/merge_default.fq - <<< "$gz_out" > /dev/null 2>&1; then
    _log_pass "gzip output matches plain reference after gunzip"
else
    _log_fail "gzip output differs from plain reference"
fi

# Smoke: -h exits 0 and prints version
run_ok "merge -h exits 0" $BFQ merge -h
run_grep "merge -h prints version" "bfqutils v" $BFQ merge -h

# Error: no input files
run_fail "merge with no input exits non-zero" $BFQ merge

# Error: mismatched R1/R2 file counts
run_fail "merge mismatched file counts" \
    $BFQ merge fixtures/testR1.fq.gz,fixtures/testR1.fq.gz fixtures/testR2.fq.gz

summary
