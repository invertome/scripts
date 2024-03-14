#!/bin/bash

# Check if a transcriptID was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <transcriptID>"
    exit 1
fi

# The string you're looking for, provided as the first command line argument
transcriptID="$1"

# Output TSV file will be named according to transcriptID
output_file=${transcriptID}.tsv

# Write the header to the output file with a tab separator
printf "Source Path\tTranscriptID\tLength\tEffectiveLength\tTPM\tNumReads\n" > "$output_file"

# Find all 'quant.sf' files in subdirectories
find . -type f -name 'quant.sf' | while read file; do
    # Search for transcriptID in each file and if found, append the line to the output file
    grep "$transcriptID" "$file" | while read matching_line; do
        # Use printf to format the output with tabs
        printf "%s\t%s\n" "$file" "$matching_line" >> "$output_file"
    done
done

echo "Done. Matched transcriptIDs are written in '$output_file'."
