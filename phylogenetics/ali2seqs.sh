#!/bin/bash

while getopts ":i:o:" opt; do
  case $opt in
    i) input_file="$OPTARG"
    ;;
    o) output_file="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

if [[ ! -f "$input_file" ]]; then
  echo "Input file not found"
  exit 1
fi

awk '/^>/ { if (NR>1) { printf("\n%s\n", $0); } else { print $0; } } /^[^>]/ { gsub("-", ""); printf("%s", $0); } END { printf("\n"); }' "$input_file" > "$output_file"
