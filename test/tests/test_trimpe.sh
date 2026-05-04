#!/bin/bash
# Tests for bfqutils trimpe
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source lib.sh
BFQ=$BFQUTILS

echo "=== trimpe ==="

# Regression: default R1 output unchanged
run_diff "default trimpe R1" \
    expected/trimpe_default_R1.fq \
    $BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
        fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Regression: default R2 output unchanged
run_diff "default trimpe R2" \
    expected/trimpe_default_R2.fq \
    $BFQ trimpe -q -1 /dev/null -2 /dev/stdout \
        fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Overlap-based trimming shortens reads (testR1/R2 have short inserts)
trim_len=$($BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
    fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null \
    | awk 'NR%4==2{sum+=length($0)} END{print sum}')
orig_len=$(gunzip -c fixtures/testR1.fq.gz | awk 'NR%4==2{sum+=length($0)} END{print sum}')
if [[ "$trim_len" -lt "$orig_len" ]]; then
    _log_pass "overlap trimming shortens R1 reads"
else
    _log_fail "overlap trimming: expected shorter reads, got same length"
fi

# No-overlap pair passes through at full input length
noop_len=$($BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
    fixtures/trimpe_r1.fq fixtures/trimpe_r2.fq 2>/dev/null \
    | awk 'NR==2{print length($0)}')
if [[ "$noop_len" -eq 20 ]]; then
    _log_pass "non-overlapping read pair preserved at input length"
else
    _log_fail "non-overlapping pair: expected len 20, got $noop_len"
fi

# Drop-both default: read2's R2 is 5bp (<15 minLen), so both mates dropped
drop_n=$($BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
    fixtures/trimpe_r1.fq fixtures/trimpe_r2.fq 2>/dev/null | grep -c '^@' || true)
if [[ "$drop_n" -eq 1 ]]; then
    _log_pass "drop-both: failed-pair mates dropped (1/2 pairs kept)"
else
    _log_fail "drop-both: expected 1 read in R1 output, got $drop_n"
fi

# Singletons via -s: surviving R1 of failed pair written to singletons file
sing_n=$($BFQ trimpe -q -1 /dev/null -2 /dev/null -s /dev/stdout \
    fixtures/trimpe_r1.fq fixtures/trimpe_r2.fq 2>/dev/null | grep -c '^@' || true)
if [[ "$sing_n" -eq 1 ]]; then
    _log_pass "singletons: surviving R1 of failed pair written to -s file"
else
    _log_fail "singletons: expected 1 singleton, got $sing_n"
fi

# Interleaved mode: line count = R1 + R2, records alternate R1/R2
r1_lines=$($BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
    fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | wc -l)
r2_lines=$($BFQ trimpe -q -1 /dev/null -2 /dev/stdout \
    fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | wc -l)
il_lines=$($BFQ trimpe -q -i \
    fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | wc -l)
if [[ "$il_lines" -eq $((r1_lines + r2_lines)) ]]; then
    _log_pass "interleaved output line count = R1 + R2"
else
    _log_fail "interleaved: expected $((r1_lines + r2_lines)) lines, got $il_lines"
fi

# Comma-separated inputs: two pairs -> output doubles
single_r1=$($BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
    fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | wc -l)
double_r1=$($BFQ trimpe -q -1 /dev/stdout -2 /dev/null \
    fixtures/testR1.fq.gz,fixtures/testR1.fq.gz \
    fixtures/testR2.fq.gz,fixtures/testR2.fq.gz 2>/dev/null | wc -l)
if [[ $((single_r1 * 2)) -eq "$double_r1" ]]; then
    _log_pass "comma-separated inputs double output"
else
    _log_fail "comma-separated: expected $((single_r1 * 2)) lines, got $double_r1"
fi

# Gzip output (-z): decompress and diff against plain reference
gz_r1=$($BFQ trimpe -q -z -1 /dev/stdout -2 /dev/null \
    fixtures/testR1.fq.gz fixtures/testR2.fq.gz 2>/dev/null | gunzip)
if diff -q expected/trimpe_default_R1.fq - <<< "$gz_r1" > /dev/null 2>&1; then
    _log_pass "gzip R1 output matches plain reference after gunzip"
else
    _log_fail "gzip R1 output differs from plain reference"
fi

# Error: neither -1/-2 nor -i given
run_fail "trimpe with no output flag exits non-zero" \
    $BFQ trimpe fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Error: -i together with -1
run_fail "trimpe -i with -1 exits non-zero" \
    $BFQ trimpe -i -1 /dev/null \
        fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Error: mismatched R1/R2 file counts
run_fail "trimpe mismatched file counts" \
    $BFQ trimpe -i \
        fixtures/testR1.fq.gz,fixtures/testR1.fq.gz fixtures/testR2.fq.gz

# Smoke: -h exits 0 and prints version
run_ok "trimpe -h exits 0" $BFQ trimpe -h
run_grep "trimpe -h prints version" "bfqutils v" $BFQ trimpe -h

# Error: no input files
run_fail "trimpe with no input exits non-zero" $BFQ trimpe

summary
