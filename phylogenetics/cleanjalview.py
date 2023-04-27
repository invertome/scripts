import argparse
from Bio import SeqIO

def clean_fasta_header(header):
    if "/" in header:
        slash_index = header.index("/")
        space_index = header.find(" ", slash_index)
        if space_index == -1:  # If there's no space after the "/"
            cleaned_header = header[:slash_index]
        else:
            cleaned_header = header[:slash_index] + header[space_index:]
        return cleaned_header
    return header

def main(input_fasta, output_fasta):
    with open(input_fasta, "r") as fasta_file, open(output_fasta, "w") as output_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            cleaned_header = clean_fasta_header(record.description)
            record.description = cleaned_header
            record.id = cleaned_header.split()[0]
            SeqIO.write(record, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean sequence headers in a FASTA file and output a new FASTA file.")
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output FASTA file", required=True)

    args = parser.parse_args()
    main(args.input, args.output)
