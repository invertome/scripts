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
parser.add_argument('-f', '--format', choices=['meme', 'pfm', 'hmm', 'psiblast', 'all'], default='meme', help='Output format (default: meme)')
parser.add_argument('-n', '--name', default='unknown', help='Motif name (default: unknown)')
parser.add_argument('-t', '--type', choices=['nt', 'aa'], required=True, help='Sequence type (nt or aa)')
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
    file_base = os.path.splitext(os.path.basename(input_file))[0]
    if args.format == 'all':
        formats = ['meme', 'pfm', 'hmm', 'psiblast']
    else:
        formats = [args.format]

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

    # Create a temporary file to hold the input alignment in FASTA format
    with NamedTemporaryFile(mode='w', delete=False, dir=args.output) as tmp_file:
        AlignIO.write(alignment, tmp_file, 'fasta')
        tmp_file.flush()

        # Run fasta-get-markov to compute background frequencies using the Docker container
        with NamedTemporaryFile(mode='w', delete=False, dir=args.output) as bg_freq_file:
            subprocess.run(['docker', 'run', '-v', f'{os.path.abspath(args.output)}:/home/meme/output', 'memesuite/memesuite', 'fasta-get-markov', os.path.join('/home/meme/output', os.path.basename(tmp_file.name)), os.path.join('/home/meme/output', os.path.basename(bg_freq_file.name))], check=True)
            bg_freq_file.flush()

            for fmt in formats:
                output_file = os.path.join(args.output, file_base + "." + fmt)

                if fmt == 'meme':
                    # Run MEME to create a .meme file using the Docker container
                    subprocess.run(['docker', 'run', '-v', f'{os.path.abspath(args.output)}:/home/meme/output', 'memesuite/memesuite', 'meme', os.path.join('/home/meme/output', os.path.basename(tmp_file.name)), '-oc', os.path.join('/home/meme/output', os.path.basename(output_file)), '-nostatus', '-bfile', os.path.join('/home/meme/output', os.path.basename(bg_freq_file.name)), '-dna' if args.type == 'nt' else '-protein'], check=True)

                elif fmt == 'pfm':
                    # Compute the frequency matrix
                    freq_matrix = count_per_position(alignment, alphabet)
                    num_seqs = len(alignment)
                    seq_length = alignment.get_alignment_length()

                    # Write the PFM file header and frequency matrix
                    with open(output_file, 'w') as outfile:
                        outfile.write("# PFM file\n")
                        outfile.write("\n")
                        outfile.write("# {}\n".format(args.name))
                        outfile.write("\n")
                        for i, symbol in enumerate(alphabet):
                            outfile.write("{}\t{}\n".format(symbol, "\t".join(["{:.2f}".format(c) for c in [row[i] for row in freq_matrix]])))
                        outfile.write("\n")

                elif fmt == 'hmm':
                    # Create a temporary file to hold the input alignment in Stockholm format
                    with NamedTemporaryFile(mode='w', delete=False, dir=args.output) as tmp_stockholm_file:
                        AlignIO.write(alignment, tmp_stockholm_file, 'stockholm')
                        tmp_stockholm_file.flush()

                        # Run hmmbuild to generate the HMM profile
                        subprocess.run(['hmmbuild', '-n', args.name, output_file, tmp_stockholm_file.name], check=True)

                elif fmt == 'psiblast':
                    # Run PSI-BLAST to create a checkpoint file
                    subprocess.run(['psiblast', '-subject', tmp_file.name, '-in_msa', tmp_file.name, '-out_ascii_pssm', output_file], check=True)

                print("Wrote {}.".format(output_file))

        # Remove temporary files
        os.remove(tmp_file.name)
        os.remove(bg_freq_file.name)
