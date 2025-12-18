#!/bin/bash
#
# Usage: ./mpf_from_string.sh <number_string> <output_file> <num_coeffs>
#
# This script finds the Julia executable, runs hpf_converter.jl to parse
# the input string, and writes the resulting coefficients and exponent
# to the specified output file.

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <number_string> <output_file> <num_coeffs>" >&2
    exit 1
fi

INPUT_STRING="$1"
OUTPUT_FILE="$2"
NUM_COEFFS="$3"

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
"$JULIA_EXE" -e "include(\"$SCRIPT_DIR/mpf_converter.jl\"); s=\"$INPUT_STRING\"; nc=$NUM_COEFFS; m,e = parse_mpf_string(s, nc); sign_int = m.is_negative ? 1 : 0; open(\"$OUTPUT_FILE\", \"w\") do f; println(f, join(vcat(m.coeffs, e, sign_int), \",\")); end" > "$OUTPUT_FILE" 2>&1