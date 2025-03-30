#!/bin/bash

echo "Running test of bfqmerge output."

if ! ../bfqmerge -q testR1.fq.gz testR2.fq.gz > test.fq ; then
    echo "Test failed, bfqmerge encountered an error."
    exit 1
fi

diff merged.fq test.fq > diff.txt

if [ -s diff.txt ] ; then
    echo "Test failed, found the following diff in bfqmerge output:"
    cat diff.txt
    exit 1
else
    echo "Test succeeded, no changes in bfqmerge output."
    rm -f test.fq
    rm -f diff.txt
    exit 0
fi

