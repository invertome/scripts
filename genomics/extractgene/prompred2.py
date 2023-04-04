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


def parse_fasta(file_path):
    """
    Parse a FASTA file and return a list of sequences.

    :param file_path: Path to the FASTA file
    :return: List of SeqRecords
    """
    with open(file_path, 'r') as fasta_file:
        return list(SeqIO.parse(fasta_file, 'fasta'))


def predict_promoter_regions(sequence, motif_file, threshold, output_path):
    """
    Use FIMO in a Docker container to predict promoter regions in the given sequence based on the provided motif file.

    :param sequence: A SeqRecord object containing the sequence to analyze
    :param motif_file: Path to the motif file
    :param threshold: FIMO threshold value
    :param output_path: Path to the output directory
    :return: List of (start, end) tuples representing promoter regions
    """
    # Create temporary directory inside the output folder
    tempdir = tempfile.mkdtemp(dir=output_path)

    # Write the sequence to a temporary FASTA file
    with open(os.path.join(tempdir, "temp_sequence.fasta"), "w") as temp_file:
        SeqIO.write(sequence, temp_file, "fasta")

    # Run FIMO within the Docker container
    fimo_output_path = os.path.join(output_path, "fimo_out")
    subprocess.run([
        "docker", "run", "-v", f"{os.getcwd()}:/home/meme",
        "memesuite/memesuite",
        "fimo",
        "--thresh", str(threshold),
        "--oc", f"/home/meme/{fimo_output_path}",
        "/home/meme/" + motif_file,
        os.path.join(tempdir, "temp_sequence.fasta")
    ])

    os.remove(os.path.join(tempdir, "temp_sequence.fasta"))

    promoter_regions = []

    with open(os.path.join(fimo_output_path, "fimo.tsv"), "r") as fimo_output_file:
        fimo_output_file.readline()  # Skip header
        for line in fimo_output_file:
            columns = line.strip().split("\t")
            if len(columns) >= 5:
                start, end = int(columns[3]) - 1, int(columns[4])
                promoter_regions.append((start, end))

    # Remove the temporary directory
    shutil.rmtree(tempdir)

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


def draw_sequence_graphics(sequences, promoter_regions_list, output_path, output_pdf, no_motif_file, output_dir):
    fig, axes = plt.subplots(len(sequences), 1, figsize=(10, len(sequences) * 2))
    if len(sequences) == 1:
        axes = [axes]

    with open(no_motif_file, "w") as no_motif:
        for ax, sequence, promoter_regions in zip(axes, sequences, promoter_regions_list):
            if promoter_regions:
                for i, (start, end) in enumerate(promoter_regions):
                    motif_sequence = str(sequence.seq[start:end])
                    label_text = f"Motif {i+1}: {start}-{end}, {motif_sequence}"
                    ax.axvspan(start, end, color=plt.cm.viridis(float(i) / len(promoter_regions)), alpha=0.5, label=label_text)
                ax.legend(loc='upper right', fontsize='small')
            else:
                no_motif.write(f"{sequence.id}\n")

            ax.set_xlim(0, len(sequence))
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_title(sequence.id)

        plt.tight_layout()
        plt.savefig(output_path)
        
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig)
        
        plt.close()


    # Write the promoter sequences to the output FASTA file
    os.makedirs(output_dir, exist_ok=True)
    output_fasta = os.path.join(output_dir, "motif_sequences.fasta")
    write_promoter_fasta(promoter_sequences, output_fasta)

    # Draw sequence graphics
    output_graphics = os.path.join(output_dir, "output_graphics.png")
    output_pdf = os.path.join(output_dir, "output_graphics.pdf")
    no_motif_file = os.path.join(output_dir, "seqs_no_motifs.list")
    draw_sequence_graphics(sequences, all_promoter_regions, output_graphics, output_pdf, no_motif_file, output_dir)

    print("Done!")
