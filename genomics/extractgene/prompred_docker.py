import os
import sys
import argparse
import subprocess
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm

def parse_fasta(file_path):
    """
    Parse a FASTA file and return a list of sequences.

    :param file_path: Path to the FASTA file
    :return: List of SeqRecords
    """
    with open(file_path, 'r') as fasta_file:
        return list(SeqIO.parse(fasta_file, 'fasta'))

def predict_promoter_regions(sequence, motif_file, threshold):
    """
    Use FIMO in a Docker container to predict promoter regions in the given sequence based on the provided motif file.

    :param sequence: A SeqRecord object containing the sequence to analyze
    :param motif_file: Path to the motif file
    :return: List of (start, end) tuples representing promoter regions
    """
    # Write the sequence to a temporary FASTA file
    with open("temp_sequence.fasta", "w") as temp_file:
        SeqIO.write(sequence, temp_file, "fasta")

    # Run FIMO within the Docker container
    subprocess.run([
        "docker", "run", "-v", f"{os.getcwd()}:/home/meme",
        "memesuite/memesuite",
        "fimo",
        "--thresh", str(threshold),
        "--oc", "/home/meme/fimo_out",
        "/home/meme/" + motif_file,
        "/home/meme/temp_sequence.fasta"
    ])

    os.remove("temp_sequence.fasta")


    # Parse the FIMO output and extract the promoter regions
    promoter_regions = []
    with open("fimo_out/fimo.tsv", "r") as fimo_output_file:
        for line in fimo_output_file:
            # Skip comments and headers
            if line.startswith("#") or line.startswith("pattern") or line.startswith("motif_id"):
                continue

            values = line.strip().split("\t")
            if len(values) >= 5:
                _, _, _, start, end, *_ = values
                promoter_regions.append((int(start), int(end)))


    return promoter_regions

def write_promoter_fasta(sequences, output_path):
    """
    Write a FASTA file containing the promoter sequences.

    :param sequences: List of SeqRecords containing promoter sequences
    :param output_path: Path to the output FASTA file
    """
    with open(output_path, "w") as output_file:
        for seq in sequences:
            SeqIO.write(seq, output_file, "fasta")
            
def draw_sequence_graphics(sequences, promoter_regions, output_path):
    gd_diagram = GenomeDiagram.Diagram("Promoter regions & Motifs")
    max_len = max([len(seq) for seq in sequences])

    for i, (seq, regions) in enumerate(zip(sequences, promoter_regions)):
        gd_track = gd_diagram.new_track(i + 1, greytrack=True, start=0, end=len(seq), scale_ticks=1, scale_smalltick_interval=1000, name=seq.id)
        gd_feature_set = gd_track.new_set()

        for start, end in regions:
            feature = SeqFeature(FeatureLocation(start, end), strand=1)
            gd_feature_set.add_feature(feature, color=colors.red, label=True, label_position="start", label_size=8)

    gd_diagram.draw(format="linear", pagesize="A4", fragments=1, start=0, end=max_len)
    gd_diagram.write(output_path, "pdf")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output directory", required=True)
    parser.add_argument("-m", "--motif", help="Path to the motif file", required=True)
    parser.add_argument("-t", "--threshold", help="FIMO threshold value", required=True, type=float)
    args = parser.parse_args()

    # Load sequences from the input FASTA file
    sequences = parse_fasta(args.input)

    # Predict and extract promoter regions for each sequence
    promoter_sequences = []
    all_promoter_regions = []  # Added to store promoter regions for drawing graphics
    for seq in sequences:
        promoter_regions = predict_promoter_regions(seq, args.motif, args.threshold)
        all_promoter_regions.append(promoter_regions)  # Store promoter regions
        for start, end in promoter_regions:
            promoter_seq = SeqRecord(seq.seq[start:end], id=f"{seq.id}_motif_{start}-{end}", description=f"Motif in promoter region for {seq.id} at position {start}-{end}")
            promoter_sequences.append(promoter_seq)

    # Write the promoter sequences to the output FASTA file
    output_fasta = os.path.join(args.output, "motif_sequences.fasta")
    os.makedirs(args.output, exist_ok=True)
    write_promoter_fasta(promoter_sequences, output_fasta)
    
    # Draw sequence graphics and save as a PDF file
    output_graphics = os.path.join(args.output, "promoter_regions_graphics.pdf")
    draw_sequence_graphics(sequences, all_promoter_regions, output_graphics)
