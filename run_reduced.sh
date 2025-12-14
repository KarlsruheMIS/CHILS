#!/bin/bash

set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_graph_file>"
    exit 1
fi

INPUT_FILE=$1
BASENAME=$(basename "$INPUT_FILE" .graph)

# If the input file path is relative, make it absolute.
if [[ "$INPUT_FILE" != /* ]]; then
    INPUT_FILE="$(pwd)/$INPUT_FILE"
fi

# Make sure we are in the script's directory, so we can find the executables.
cd "$(dirname "$0")"

LR_OUT=$(./MWIS_REDUCE "$INPUT_FILE" --time_limit=1800 --kernel="$BASENAME" | tr -d '\n')

rm "$BASENAME"_weight.csv
rm "$BASENAME".csv

LR_TIME=$(echo "$LR_OUT" | awk -F',' '{print $11}')

CHILS_TIMEOUT=$(echo "3600 - $LR_TIME" | bc)

CHILS_OUT=$(./MWIS_CHILS -g "$INPUT_FILE" -t "$CHILS_TIMEOUT" -c 1)

rm "$BASENAME".kernel_graph

echo "$LR_OUT,$CHILS_OUT"