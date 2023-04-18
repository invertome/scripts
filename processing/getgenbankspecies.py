import argparse
from Bio import SeqIO

def parse_fasta_header(header):
    if "[" in header and "]" in header:
        identifier = header.split()[0].strip(">")
        species_start = header.index("[") + 1
        species_end = header.index("]")
        species = header[species_start:species_end].replace(" ", "_")
        return identifier, species
    return None

def main(input_fasta, output_tsv):
    with open(input_fasta, "r") as fasta_file, open(output_tsv, "w") as tsv_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            print(f"Processing header: {record.description}")  # Debug statement
            parsed_header = parse_fasta_header(record.description)
            if parsed_header:
                identifier, species = parsed_header
                tsv_file.write(f"{identifier}\t{species}\n")
            else:
                print(f"Skipped header: {record.description}")  # Debug statement

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a FASTA file and output a TSV file with specific fields.")
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output TSV file", required=True)

    args = parser.parse_args()
    main(args.input, args.output)
