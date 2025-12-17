#!/bin/bash
#
# Usage: ./mpf_to_string.sh <c0> <c1> <c2> <c3> <exp> <sign> <output_file>
#
# This script finds the Julia executable, runs hpf_converter.jl to convert
# the given components into a string, and writes the result to the output file.

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <c0> <c1> <c2> <c3> <exp> <sign> <output_file>" >&2
    exit 1
fi

C0="$1"; C1="$2"; C2="$3"; C3="$4"; EXP="$5"; SIGN="$6"; OUTPUT_FILE="$7"

# --- Find Julia Executable (same logic as check_julia_version.sh) ---
JULIA_EXE=""
if [ -n "$JULIA_PATH" ] && [ -x "$JULIA_PATH" ]; then
    JULIA_EXE="$JULIA_PATH"
elif command -v julia >/dev/null 2>&1; then
    JULIA_EXE=$(command -v julia)
fi
if [ -z "$JULIA_EXE" ] && [ -d "$HOME/julia" ]; then
    LATEST_JULIA_DIR=$(ls -v "$HOME/julia" | grep '^julia-' | tail -n 1)
    if [ -n "$LATEST_JULIA_DIR" ]; then
        CANDIDATE_PATH="$HOME/julia/$LATEST_JULIA_DIR/bin/julia"
        if [ -x "$CANDIDATE_PATH" ]; then
            JULIA_EXE="$CANDIDATE_PATH"
        fi
    fi
fi

if [ -z "$JULIA_EXE" ]; then
    echo "ERROR: Julia executable not found." > "$OUTPUT_FILE"
    exit 1
fi

# Construct and execute the Julia command.
SCRIPT_DIR=$(dirname -- "$0")
JULIA_CMD="include(\"$SCRIPT_DIR/mpf_converter.jl\"); coeffs=Int64[$C0,$C1,$C2,$C3]; e=$EXP; is_neg=$SIGN==1; m=MPF(coeffs,Int32(e),is_neg); result=to_string(m); open(\"$OUTPUT_FILE\",\"w\") do f; println(f,result); end"

"$JULIA_EXE" -e "$JULIA_CMD" > /dev/null 2>&1