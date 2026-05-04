#!/bin/bash
# Test helper library for bfqutils.
# Run test.sh from the test/ directory; all paths are relative to test/.
#
# To regenerate reference outputs after an intentional behaviour change:
#   (cd test && bash lib.sh regen)
# Each test_*.sh must define a regen() function and check $REGEN at the top.

BFQUTILS=../bfqutils
PASS=0
FAIL=0
REGEN=${REGEN:-0}

_log_pass() { echo "  PASS  $1"; ((PASS++)); }
_log_fail() { echo "  FAIL  $1"; ((FAIL++)); }

# run_diff <name> <expected_file> <command...>
# Runs command, diffs stdout against expected_file.
run_diff() {
    local name=$1 expected=$2; shift 2
    local out
    if out=$("$@" 2>/dev/null); then
        if diff -q "$expected" - <<< "$out" > /dev/null 2>&1; then
            _log_pass "$name"
        else
            _log_fail "$name"
            diff "$expected" - <<< "$out" | head -20
        fi
    else
        _log_fail "$name (command failed with exit $?)"
    fi
}

# run_diff_stderr <name> <expected_file> <command...>
# Runs command, diffs stderr against expected_file.
run_diff_stderr() {
    local name=$1 expected=$2; shift 2
    local err
    if err=$("$@" 2>&1 >/dev/null); then
        if diff -q "$expected" - <<< "$err" > /dev/null 2>&1; then
            _log_pass "$name"
        else
            _log_fail "$name"
            diff "$expected" - <<< "$err" | head -20
        fi
    else
        _log_fail "$name (command failed with exit $?)"
    fi
}

# run_ok <name> <command...>
# Passes if command exits 0.
run_ok() {
    local name=$1; shift
    if "$@" > /dev/null 2>&1; then
        _log_pass "$name"
    else
        _log_fail "$name (expected exit 0, got $?)"
    fi
}

# run_fail <name> <command...>
# Passes if command exits non-zero.
run_fail() {
    local name=$1; shift
    if "$@" > /dev/null 2>&1; then
        _log_fail "$name (expected non-zero exit)"
    else
        _log_pass "$name"
    fi
}

# run_grep <name> <pattern> <command...>
# Passes if stdout+stderr contains pattern.
run_grep() {
    local name=$1 pattern=$2; shift 2
    local out
    out=$("$@" 2>&1)
    if echo "$out" | grep -q "$pattern"; then
        _log_pass "$name"
    else
        _log_fail "$name (pattern '$pattern' not found in output)"
        echo "$out" | head -5
    fi
}

# run_no_match <name> <pattern> <command...>
# Passes if stdout does NOT contain pattern.
run_no_match() {
    local name=$1 pattern=$2; shift 2
    local out
    out=$("$@" 2>/dev/null)
    if echo "$out" | grep -q "$pattern"; then
        _log_fail "$name (pattern '$pattern' found in output)"
        echo "$out" | grep "$pattern" | head -3
    else
        _log_pass "$name"
    fi
}

summary() {
    local total=$((PASS + FAIL))
    echo ""
    echo "Results: $PASS/$total passed"
    [[ $FAIL -eq 0 ]]
}
