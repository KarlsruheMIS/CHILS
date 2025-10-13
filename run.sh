#!/bin/bash

set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_graph_file>"
    exit 1
fi

INPUT_FILE=$1
BASENAME=$(basename "$INPUT_FILE" .graph)
JSON_FILE="${BASENAME}.json"
OUTPUT_DIR="${BASENAME}"

# Make sure we are in the script's directory, so we can find the executables.
cd "$(dirname "$0")"

# If the input file path is relative, make it absolute.
if [[ "$INPUT_FILE" != /* ]]; then
    INPUT_FILE="$(pwd)/$INPUT_FILE"
fi

# echo "Processing $INPUT_FILE..."

# echo "1. Converting to JSON..."
METIS_OUT=$(./METIS_TO_JSON "$INPUT_FILE" "$JSON_FILE")

# echo "2. Creating output directory..."
mkdir -p "$OUTPUT_DIR"

# echo "3. Running BSA_GREEDY..."
BSA_OUT=$(./BSA_GREEDY "$JSON_FILE" 16)

METIS_TIME=$(echo "$METIS_OUT" | awk -F',' '{print $3}')
BSA_TIME=$(echo "$BSA_OUT" | awk -F',' '{print $1}')
CHILS_TIMEOUT=$(echo "3600 - $METIS_TIME - $BSA_TIME" | bc)

# echo "4. Running CHILS for $CHILS_TIMEOUT seconds..."
CHILS_OUT=$(./MWIS_CHILS -g "$INPUT_FILE" -f "$OUTPUT_DIR" -t "$CHILS_TIMEOUT" -c 1)

# echo "5. Cleaning up..."
rm -rf "$JSON_FILE"
rm -rf "$OUTPUT_DIR"

# echo "Done."

echo "$METIS_OUT,$BSA_OUT,$CHILS_OUT"