#!/bin/bash
# Tests for bfqutils stats
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source lib.sh
BFQ=$BFQUTILS

echo "=== stats ==="

# Regression: default summary output (stderr) unchanged
run_diff_stderr "default stats summary" \
    expected/stats_default.txt \
    $BFQ stats fixtures/testR1.fq.gz

# All histogram outputs: each diff against a committed reference
run_diff "stats -l length histogram" \
    expected/stats_lenhist.tsv \
    $BFQ stats -N -l /dev/stdout fixtures/testR1.fq.gz

run_diff "stats -g GC histogram" \
    expected/stats_gchist.tsv \
    $BFQ stats -N -g /dev/stdout fixtures/testR1.fq.gz

run_diff "stats -q mean quality histogram" \
    expected/stats_qualhist.tsv \
    $BFQ stats -N -q /dev/stdout fixtures/testR1.fq.gz

run_diff "stats -Q per-position quality" \
    expected/stats_ppqual.tsv \
    $BFQ stats -N -Q /dev/stdout fixtures/testR1.fq.gz

run_diff "stats -b per-position base content" \
    expected/stats_ppbase.tsv \
    $BFQ stats -N -b /dev/stdout fixtures/testR1.fq.gz

run_diff "stats -k kmer counts" \
    expected/stats_kmers.tsv \
    $BFQ stats -N -k /dev/stdout fixtures/testR1.fq.gz

run_diff "stats -K 4 kmer counts" \
    expected/stats_kmers_K4.tsv \
    $BFQ stats -N -K 4 -k /dev/stdout fixtures/testR1.fq.gz

# -n 50: only first 50 reads processed
n50_reads=$($BFQ stats -n 50 fixtures/testR1.fq.gz 2>&1 | grep 'Total reads' | grep -o '[0-9]*')
if [[ "$n50_reads" -eq 50 ]]; then
    _log_pass "-n 50 processes exactly 50 reads"
else
    _log_fail "-n 50: expected 50 reads, stats reported $n50_reads"
fi

# Bug regression (v1.2.1): Median must reflect read count, not base count.
# median_known.fq has 5×50bp + 5×100bp reads → median should be 50.
median=$($BFQ stats fixtures/median_known.fq 2>&1 | grep 'Median length' | grep -o '[0-9]*')
if [[ "$median" -eq 50 ]]; then
    _log_pass "median length is 50 for 5×50bp + 5×100bp input"
else
    _log_fail "median: expected 50, got $median"
fi

# -O passthrough: output bytes equal input after decompression
in_bytes=$(gunzip -c fixtures/testR1.fq.gz | wc -c)
out_bytes=$($BFQ stats -N -O fixtures/testR1.fq.gz 2>/dev/null | wc -c)
if [[ "$in_bytes" -eq "$out_bytes" ]]; then
    _log_pass "-O passthrough byte count matches input"
else
    _log_fail "-O passthrough: input $in_bytes bytes, output $out_bytes bytes"
fi

# Bug regression (v1.2.3): -O on comment-less headers must not print "(null)".
run_no_match "-O no_comment.fq has no (null)" "(null)" \
    $BFQ stats -N -O fixtures/no_comment.fq

# Bug regression (v1.2.2): empty reads must not produce NaN or UB.
run_ok "empty_reads.fq exits cleanly" \
    $BFQ stats fixtures/empty_reads.fq
run_no_match "empty reads: no 'nan' in output" "nan" \
    $BFQ stats fixtures/empty_reads.fq

# Bug regression (v1.2.2): all-N reads must give obs/exp = 0 (not nan/inf).
run_ok "all_n.fq exits cleanly" \
    $BFQ stats fixtures/all_n.fq
run_no_match "all-N reads: no nan/inf in obs/exp" $'\(nan\|inf\)' \
    $BFQ stats fixtures/all_n.fq

# Bug regression (v1.2.2): invalid -n argument rejected.
run_fail "-n abc rejected" $BFQ stats -n abc fixtures/testR1.fq.gz
run_grep "-n abc prints error" "Error" $BFQ stats -n abc fixtures/testR1.fq.gz

# Bug regression (v1.2.2): -z without -O rejected.
run_fail "-z without -O exits non-zero" $BFQ stats -z fixtures/testR1.fq.gz
run_grep "-z without -O mentions -O" "\-O" $BFQ stats -z fixtures/testR1.fq.gz

# Stdin (-) input: summary should match normal file run
stdin_out=$($BFQ stats - < fixtures/testR1.fq.gz 2>&1)
file_out=$($BFQ stats fixtures/testR1.fq.gz 2>&1)
if [[ "$stdin_out" == "$file_out" ]]; then
    _log_pass "stdin (-) input matches file input"
else
    _log_fail "stdin (-) output differs from file input"
fi

# Smoke: -h exits 0 and prints version
run_ok "stats -h exits 0" $BFQ stats -h
run_grep "stats -h prints version" "bfqutils v" $BFQ stats -h

# Error: no input files
run_fail "stats with no input exits non-zero" $BFQ stats

summary
