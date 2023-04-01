import os
import sys
import argparse
import subprocess
import shutil
import tempfile
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import Diagram as _Diagram, _LinearDrawer, CrossLink
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.units import cm
from reportlab.graphics.shapes import String
from reportlab.pdfgen import canvas



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
    fimo_output_path = os.path.join(args.output, "fimo_out")
    subprocess.run([
        "docker", "run", "-v", f"{os.getcwd()}:/home/meme",
        "memesuite/memesuite",
        "fimo",
        "--thresh", str(threshold),
        "--oc", f"/home/meme/{fimo_output_path}",
        "/home/meme/" + motif_file,
        "/home/meme/temp_sequence.fasta"
    ])

    os.remove("temp_sequence.fasta")

    promoter_regions = []

    with open(os.path.join(fimo_output_path, "fimo.tsv"), "r") as fimo_output_file:
        fimo_output_file.readline()  # Skip header
        for line in fimo_output_file:
            columns = line.strip().split("\t")
            if len(columns) >= 5:
                start, end = int(columns[3]) - 1, int(columns[4])
                promoter_regions.append((start, end))

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
            
class CustomDiagram(_Diagram):
    def new_track(self, track_level, *args, **kwargs):
        track = super().new_track(track_level, *args, **kwargs)
        track.height = 1.0
        return track

def draw_sequence_graphics(sequences, promoter_regions_list, output_path):
    max_len = max(len(seq) for seq in sequences)
    pdf_canvas = canvas.Canvas(output_path, pagesize=letter)
    
    for sequence, promoter_regions in zip(sequences, promoter_regions_list):
        gd_diagram = GenomeDiagram.Diagram(sequence.id)
        feature_track = gd_diagram.new_track(1, name=f"Track_{sequence.id}", greytrack=True, greytrack_labels=10, scale=True, height=1.0, start_pad=0.5, end_pad=0.5)
        feature_set = feature_track.new_set()

        for i, (start, end) in enumerate(promoter_regions):
            feature = SeqFeature(FeatureLocation(start, end), strand=1)
            feature_set.add_feature(
                feature,
                color=colors.Color(*plt.cm.viridis(float(i) / len(promoter_regions))[:3]),
                name=f"Motif {i+1}",
                label=True,
                label_size=8,
                label_color=colors.black,
                label_position="middle",
            )
        
        gd_diagram.draw(
            format="linear",
            orientation="landscape",
            pagesize=(max_len * 0.01 * cm, 10 * cm),
            fragments=1,
            start=0,
            end=len(sequence),
            tracklines=0,
            x=0.05,
            y=0.85,
        )
        gd_diagram.write_to_pdf(pdf_canvas)
        pdf_canvas.showPage()
    
    pdf_canvas.save()



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
        promoter_regions = predict_promoter_regions(seq, args.motif, args.threshold)
        all_promoter_regions.append(promoter_regions)  # Store promoter regions
        for start, end in promoter_regions:
            promoter_seq = SeqRecord(seq.seq[start:end], id=f"{seq.id}_motif_{start}-{end}", description=f"Motif in promoter region for {seq.id} at position {start}-{end}")
            promoter_sequences.append(promoter_seq)

    # Write the promoter sequences to the output FASTA file
    output_fasta = os.path.join(args.output, "motif_sequences.fasta")
    os.makedirs(args.output, exist_ok=True)
    write_promoter_fasta(promoter_sequences, output_fasta)
    
    # Draw sequence graphics
    draw_sequence_graphics(sequences, all_promoter_regions, output_graphics)
