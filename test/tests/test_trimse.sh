#!/bin/bash
# Tests for bfqutils trimse
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source lib.sh
BFQ=$BFQUTILS

echo "=== trimse ==="

# Regression: default trim output unchanged
run_diff "default trimse" \
    expected/trimse_default.fq \
    $BFQ trimse -q fixtures/testR1.fq.gz

# Bug regression (v1.2.1): reads with no adapter match must be kept, not dropped.
# no_adapter.fq has 5 reads with no Illumina adapter → all 5 must appear in output.
noAdp_n=$($BFQ trimse -q fixtures/no_adapter.fq 2>/dev/null | grep -c '^@')
if [[ "$noAdp_n" -eq 5 ]]; then
    _log_pass "no-adapter reads retained (5/5)"
else
    _log_fail "no-adapter reads: expected 5 in output, got $noAdp_n"
fi

# Bug regression (v1.2.1): default -M 15 filters sub-15 bp reads.
# short.fq has 2 reads <15 bp and 2 reads ≥15 bp; only the longer 2 should survive.
short_n=$($BFQ trimse -q fixtures/short.fq 2>/dev/null | grep -c '^@')
if [[ "$short_n" -eq 2 ]]; then
    _log_pass "sub-15 bp reads filtered by default -M (2/4 kept)"
else
    _log_fail "short reads: expected 2 in output, got $short_n"
fi

# -M 50 stricter filter: short.fq has no reads ≥50 bp, so output should be empty.
M50_n=$($BFQ trimse -q -M 50 fixtures/short.fq 2>/dev/null | grep -c '^@' || true)
if [[ "$M50_n" -eq 0 ]]; then
    _log_pass "-M 50 filters all reads in short.fq"
else
    _log_fail "-M 50: expected 0 reads in output, got $M50_n"
fi

# -g flag was dead code before v1.2.2: verify it is accepted and changes output.
run_ok "-g flag accepted" \
    $BFQ trimse -q -g 5 fixtures/polyg.fq

# polyg.fq has 3 reads with ≥15-G tails. With -g 5, the polyG region should be
# trimmed, shortening each read vs the untrimmed original.
orig_len=$($BFQ trimse -q -g 999 fixtures/polyg.fq 2>/dev/null \
    | awk 'NR%4==2{sum+=length($0)} END{print sum}')
trimmed_len=$($BFQ trimse -q -g 5 fixtures/polyg.fq 2>/dev/null \
    | awk 'NR%4==2{sum+=length($0)} END{print sum}')
if [[ "$trimmed_len" -lt "$orig_len" ]]; then
    _log_pass "-g 5 shortens reads with polyG tails"
else
    _log_fail "-g 5 did not shorten reads: orig=$orig_len trimmed=$trimmed_len"
fi

# Bug regression (v1.2.2): trim3p precision.
# q_borderline.fq has a Q15 tail (should NOT be trimmed with -t 15)
# and a Q14 tail (should be trimmed with -t 15).
# With the old int-cast bug, both were trimmed (Q15 → Q10 after cast < Q15).
q15_len=$($BFQ trimse -q -t 15 -w 5 fixtures/q_borderline.fq 2>/dev/null \
    | awk 'NR==2{print length($0)}')
q14_len=$($BFQ trimse -q -t 15 -w 5 fixtures/q_borderline.fq 2>/dev/null \
    | awk 'NR==6{print length($0)}')
if [[ "$q15_len" -eq 25 ]]; then
    _log_pass "Q15 tail preserved with -t 15 (trim3p precision)"
else
    _log_fail "Q15 tail trimmed with -t 15: expected len 25, got $q15_len"
fi
if [[ "$q14_len" -lt "$q15_len" ]]; then
    _log_pass "Q14 tail trimmed more than Q15 tail with -t 15 (trim3p precision)"
else
    _log_fail "Q14 tail not trimmed with -t 15: q14_len=$q14_len q15_len=$q15_len"
fi

# Comma-separated inputs: two copies → double the output
single_n=$($BFQ trimse -q fixtures/testR1.fq.gz 2>/dev/null | wc -l)
double_n=$($BFQ trimse -q fixtures/testR1.fq.gz,fixtures/testR1.fq.gz 2>/dev/null | wc -l)
if [[ $((single_n * 2)) -eq "$double_n" ]]; then
    _log_pass "comma-separated inputs double output"
else
    _log_fail "comma-separated: expected $((single_n * 2)) lines, got $double_n"
fi

# Gzip output (-z) matches plain reference after decompression
gz_out=$($BFQ trimse -q -z fixtures/testR1.fq.gz 2>/dev/null | gunzip)
if diff -q expected/trimse_default.fq - <<< "$gz_out" > /dev/null 2>&1; then
    _log_pass "gzip output matches plain reference after gunzip"
else
    _log_fail "gzip output differs from plain reference"
fi

# Stdin (-) input
stdin_out=$(gunzip -c fixtures/testR1.fq.gz | $BFQ trimse -q - 2>/dev/null)
if diff -q expected/trimse_default.fq - <<< "$stdin_out" > /dev/null 2>&1; then
    _log_pass "stdin (-) input works"
else
    _log_fail "stdin (-) output differs from file input"
fi

# Smoke: -h exits 0 and prints version
run_ok "trimse -h exits 0" $BFQ trimse -h
run_grep "trimse -h prints version" "bfqutils v" $BFQ trimse -h

# Error: no input files
run_fail "trimse with no input exits non-zero" $BFQ trimse

summary
