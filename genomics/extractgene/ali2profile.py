import argparse
import os
import subprocess
from Bio import AlignIO
from tempfile import NamedTemporaryFile
from collections import Counter

# Function to count the number of occurrences of each nucleotide at each position in the alignment
def count_per_position(alignment):
    length = alignment.get_alignment_length()
    counts = [Counter(alignment[:, i]) for i in range(length)]
    matrix = [[counts[i][nuc] for nuc in ['A', 'C', 'G', 'T']] for i in range(length)]
    return matrix

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Convert multiple sequence alignments to motif profile files.')
parser.add_argument('-i', '--input', required=True, nargs='+', help='Path(s) to input FASTA file(s)')
parser.add_argument('-o', '--output', required=True, help='Path to output folder')
parser.add_argument('-f', '--format', choices=['meme', 'pfm', 'hmm'], default='meme', help='Output format (default: meme)')
parser.add_argument('-n', '--name', default='unknown', help='Motif name (default: unknown)')
parser.add_argument('-e', '--empfreq', action='store_true', help='Use empirical nucleotide frequencies from the alignment')
parser.add_argument('-b', '--bgfreq', type=float, nargs=4, default=[0.25, 0.25, 0.25, 0.25], help='Background nucleotide frequencies (default: 0.25 for each nucleotide)')
args = parser.parse_args()

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
    num_seqs = len(alignment)

    if args.empfreq:
        # Compute empirical nucleotide frequencies
        background_freq = []
        for i in range(alignment.get_alignment_length()):
            col = alignment[:, i]
            counts = Counter(col)
            freqs = [counts[c] / num_seqs for c in ['A', 'C', 'G', 'T']]
            background_freq.append(freqs)
        background_freq = [sum(x) / len(background_freq) for x in zip(*background_freq)]
    else:
        # Use flat background nucleotide frequencies
        background_freq = args.bgfreq

    if args.format == 'meme':
        # Compute the position frequency matrix (PFM)
        freq_matrix = count_per_position(alignment)
        seq_length = len(freq_matrix)

        # Write the MEME motif file header and frequency matrix
        with open(output_file, 'w') as outfile:
            outfile.write("MEME version 4\n")
            outfile.write("\n")
            outfile.write("ALPHABET= ACGT\n")
            outfile.write("STRANDS= + -\n")
            outfile.write("Background letter frequencies (from unknown source):\n")
            outfile.write("A {:.4f} C {:.4f} G {:.4f} T {:.4f}\n".format(*background_freq))
            outfile.write("\n")
            outfile.write("MOTIF {}\n".format(args.name))
            outfile.write("letter-probability matrix: alength= 4 w= {} nsites= {}\n".format(seq_length, num_seqs))
            for i in range(seq_length):
                outfile.write("{:.4f} {:.4f} {:.4f} {:.4f}\n".format(*freq_matrix[i]))
            outfile.write("\n")

    elif args.format == 'pfm':
        # Compute the position frequency matrix (PFM)
        freq_matrix = count_per_position(alignment)

        # Write the PFM file header and frequency matrix
        with open(output_file, 'w') as outfile:
            outfile.write("# PFM file\n")
            outfile.write("\n")
            outfile.write("# {}\n".format(args.name))
            for i in range(4):
                outfile.write("# {} [{:.2f} {:.2f} {:.2f} {:.2f}]\n".format(['A', 'C', 'G', 'T'][i], *background_freq))
            outfile.write("\n")
            outfile.write("{}\n".format("\t".join(["{}".format(i+1) for i in range(len(freq_matrix))])))
            for i in range(4):
                outfile.write("{}\t{}\n".format(['A', 'C', 'G', 'T'][i], "\t".join(["{}".format(c) for c in freq_matrix[:,i]])))
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

