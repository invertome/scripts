import os
import sys
import argparse
import subprocess
import shutil
import tempfile
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# ... (parse_fasta and predict_promoter_regions functions unchanged) ...

def write_promoter_fasta(sequences, output_path):
    """
    Write a FASTA file containing the promoter sequences.

    :param sequences: List of SeqRecords containing promoter sequences
    :param output_path: Path to the output FASTA file
    """
    with open(output_path, "w") as output_file:
        for seq in sequences:
            seq.description = seq.id  # Replacing description with motif ID
            SeqIO.write(seq, output_file, "fasta")

# ... (draw_sequence_graphics function unchanged) ...

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output directory", required=True)
    parser.add_argument("-m", "--motif", help="Path to the motif file", required=True)
    parser.add_argument("-t", "--threshold", help="FIMO threshold value", required=True, type=float)
    args = parser.parse_args()

    # Create output folders
    os.makedirs(args.output, exist_ok=True)
    output_graphics = os.path.join(args.output, "graphics")
    os.makedirs(output_graphics, exist_ok=True)

    # Load sequences from the input FASTA file
    sequences = parse_fasta(args.input)

    # Predict and extract promoter regions for each sequence
    promoter_sequences = []
    all_promoter_regions = []  # Added to store promoter regions for drawing graphics
    for seq in sequences:
        seq_id = seq.id.split(':')[0]  # Use only the part before the colon for creating folders
        output_seq_folder = os.path.join(args.output, seq_id)
        os.makedirs(output_seq_folder, exist_ok=True)
        fimo_output_path = os.path.join(output_seq_folder, "fimo_out")

        promoter_regions = predict_promoter_regions(seq, args.motif, args.threshold, fimo_output_path)
        all_promoter_regions.append(promoter_regions)  # Store promoter regions
        for start, end in promoter_regions:
            promoter_seq = SeqRecord(seq.seq[start:end], id=f"{seq.id}_motif_{start}-{end}", description=f"Motif in promoter region for {seq.id} at position {start}-{end}")
            promoter_sequences.append(promoter_seq)

        # Write the promoter sequences to the output FASTA file
        output_fasta = os.path.join(args.output, "motif_sequences.fasta")
        os.makedirs(args.output, exist_ok=True)
        write_promoter_fasta(promoter_sequences, output_fasta)

        # Draw sequence graphics
        output_graphics = os.path.join(args.output, "output_graphics.png")
        output_pdf = os.path.join(args.output, "output_graphics.pdf")
        no_motif_file = os.path.join(args.output, "seqs_no_motifs.list")
        draw_sequence_graphics(sequences, all_promoter_regions, output_graphics, output_pdf, no_motif_file)
