import argparse
from Bio import SeqIO
import re

def parse_fasta_header(header):
    pattern = r'>\S+\s.*\[(\S+\s\S+)\]'
    match = re.search(pattern, header)
    if match:
        species = match.group(1).replace(" ", "_")
        return species
    return None

def main(input_fasta, output_tsv):
    with open(input_fasta, "r") as fasta_file, open(output_tsv, "w") as tsv_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            species = parse_fasta_header(record.description)
            if species:
                tsv_file.write(f"{record.id}\t{species}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a FASTA file and output a TSV file with specific fields.")
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output TSV file", required=True)

    args = parser.parse_args()
    main(args.input, args.output)
