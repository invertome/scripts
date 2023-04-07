import argparse
import os
import subprocess
from Bio import AlignIO
from tempfile import NamedTemporaryFile
from collections import Counter

# Function to count the number of occurrences of each symbol at each position in the alignment
def count_per_position(alignment, alphabet):
    length = alignment.get_alignment_length()
    counts = [Counter(alignment[:, i]) for i in range(length)]
    matrix = [[counts[i][symbol] for symbol in alphabet] for i in range(length)]
    return matrix

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Convert multiple sequence alignments to motif profile files.')
parser.add_argument('-i', '--input', required=True, nargs='+', help='Path(s) to input FASTA file(s)')
parser.add_argument('-o', '--output', required=True, help='Path to output folder')
parser.add_argument('-f', '--format', choices=['meme', 'pfm', 'hmm'], default='meme', help='Output format (default: meme)')
parser.add_argument('-n', '--name', default='unknown', help='Motif name (default: unknown)')
parser.add_argument('-t', '--type', choices=['nt', 'aa'], required=True, help='Sequence type (nt or aa)')
parser.add_argument('-e', '--empfreq', action='store_true', help='Use empirical symbol frequencies from the alignment')
parser.add_argument('-b', '--bgfreq', type=float, nargs='+', default=None, help='Background symbol frequencies (default: flat frequency)')
args = parser.parse_args()

# Define the nucleotide and amino acid alphabets
nucleotides = set("ACGTUacgtu")
amino_acids = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

# Create the output folder if it doesn't exist
if not os.path.exists(args.output):
    os.makedirs(args.output)

# Process each input file
for input_file in args.input:

    # Determine the output file path based on the input file name and output format
    if args.format == 'meme':
        output_file = os.path.join(args.output, os.path.splitext(os.path.basename(input_file))[0] + ".meme")
    elif args.format == 'pfm':
        output_file = os.path.join(args.output, os.path.splitext(os.path.basename(input_file))[0] + ".pfm")
    elif args.format == 'hmm':
        output_file = os.path.join(args.output, os.path.splitext(os.path.basename(input_file))[0] + ".hmm")

    # Load the input alignment file
    alignment = AlignIO.read(input_file, 'fasta')

    # Convert U to T if it exists in the sequence for nucleotide sequences
    if args.type == 'nt':
        for seq in alignment:
            seq.seq = seq.seq.replace('U', 'T')

    # Determine the alphabet
    if args.type == 'nt':
        alphabet = ['A', 'C', 'G', 'T']
    elif args.type == 'aa':
        alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Determine the background frequency
    if args.bgfreq is not None:
        background_freq = args.bgfreq
    else:
        background_freq = [1.0 / len(alphabet)] * len(alphabet)

    # Convert U to T if it exists in the sequence
    if args.type == 'nt':
        for seq in alignment:
            seq.seq = seq.seq.replace('U', 'T')

    # Compute the frequency matrix
    if args.empfreq:
        freq_matrix = count_per_position(alignment, alphabet)
        num_seqs = len(alignment)
        seq_length = alignment.get_alignment_length()
    else:
        counts = alignment[:, :]
        num_seqs = len(alignment)
        seq_length = alignment.get_alignment_length()
        freq_matrix = [[(counts[:, i] == symbol).sum() for symbol in alphabet] for i in range(seq_length)]
        freq_matrix = [[(symbol + 0.25) / (num_seqs + 1) for symbol in row] for row in freq_matrix]

    # Write the output file
    if args.format == 'meme':
        # Write the MEME motif file header and frequency matrix
        with open(output_file, 'w') as outfile:
            outfile.write("MEME version 4\n")
            outfile.write("\n")
            outfile.write("ALPHABET= {}\n".format(''.join(alphabet)))
            outfile.write("STRANDS= + -\n")
            outfile.write("Background letter frequencies (from unknown source):\n")
            outfile.write("{}\n".format(' '.join('{:.4f}'.format(f) for f in background_freq)))
            outfile.write("\n")
            outfile.write("MOTIF {}\n".format(args.name))
            outfile.write("letter-probability matrix: alength= {} w= {} nsites= {}\n".format(len(alphabet), seq_length, num_seqs))
            for i in range(seq_length):
                outfile.write("{}\n".format(' '.join('{:.4f}'.format(f) for f in freq_matrix[i])))
            outfile.write("\n")

    elif args.format == 'pfm':
        # Write the PFM file header and frequency matrix
        with open(output_file, 'w') as outfile:
            outfile.write("# PFM file\n")
            outfile.write("\n")
            outfile.write("# {}\n".format(args.name))
            for i, symbol in enumerate(alphabet):
                outfile.write("# {} {}\n".format(symbol, ' '.join('{:.2f}'.format(f) for f in background_freq)))
            outfile.write("\n")
            outfile.write("{}\n".format("\t".join(["{}".format(i+1) for i in range(len(freq_matrix))])))
            for i, symbol in enumerate(alphabet):
                outfile.write("{}\t{}\n".format(symbol, "\t".join(["{}".format(c) for c in freq_matrix[:,i]])))
            outfile.write("\n")

    elif args.format == 'hmm':
        # Create a temporary file to hold the input alignment in Stockholm format
        with NamedTemporaryFile(mode='w', delete=False) as tmp_file:
            AlignIO.write(alignment, tmp_file, 'stockholm')
            tmp_file.flush()

            # Run hmmbuild to generate the HMM profile
            hmm_file = os.path.splitext(output_file)[0] + ".hmm"
            subprocess.run(['hmmbuild', '-n', args.name, hmm_file, tmp_file.name], check=True)

    print("Wrote {}.".format(output_file))
