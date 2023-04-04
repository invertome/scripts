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

def predict_promoter_regions(sequence, motif_file, threshold, fimo_output_path):
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

    # Copy the fimo.gff file to main output folder
    shutil.copy(os.path.join(fimo_output_path, "fimo.gff"), os.path.join(args.output, f"{sequence.id}_fimo.gff"))

    return promoter_regions

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

def draw_sequence_graphics(sequences, promoter_regions_list, output_path, output_pdf, no_motif_file):
    fig, axes = plt.subplots(len(sequences), 1, figsize=(10, len(sequences) * 2.5))  # Increase the figure height
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

            # Reduce the height of the subplot while keeping the same proportions
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width, box.height * 0.8])

    plt.subplots_adjust(hspace=0.5)  # Add spacing between subplots with the hspace parameter
    plt.savefig(output_path)

    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig)

    plt.close()



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

    fimo_output_path = os.path.join(args.output, "fimo_out")
    os.makedirs(fimo_output_path, exist_ok=True)

    # Load sequences from the input FASTA file
    sequences = parse_fasta(args.input)

    # Initialize summary_tsv and combined_gff files
    summary_tsv = os.path.join(args.output, "summary_fimo.tsv")
    combined_gff = os.path.join(args.output, "combined_fimo.gff")
    
    # Initialize summary_tsv and combined_gff files with header lines
    with open(summary_tsv, "w") as summary_file, open(combined_gff, "w") as gff_file:
        summary_file.write("#pattern_name\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n")
        gff_file.write("##gff-version 3\n")

    # Initialize variables
    all_promoter_regions = []
    promoter_sequences = []

    # Process each sequence        
    for seq in sequences:
        promoter_regions = predict_promoter_regions(seq, args.motif, args.threshold, fimo_output_path)

        # Append fimo.tsv and fimo.gff content to summary files
        with open(os.path.join(fimo_output_path, "fimo.tsv"), "r") as fimo_tsv, open(os.path.join(args.output, f"{seq.id}_fimo.gff"), "r") as fimo_gff, open(summary_tsv, "a") as summary_file, open(combined_gff, "a") as gff_file:
            fimo_tsv.readline()  # Skip header
            fimo_gff.readline()  # Skip header

            for line in fimo_tsv:
                summary_file.write(line)

            for line in fimo_gff:
                gff_file.write(line)

        all_promoter_regions.append(promoter_regions)  # Store promoter regions
        for start, end in promoter_regions:
            promoter_seq = SeqRecord(seq.seq[start:end], id=f"{seq.id}_motif_{start}-{end}", description=f"Motif in promoter region for {seq.id} at position {start}-{end}")
            promoter_sequences.append(promoter_seq)
            
    # Remove individual _fimo.gff files after creating the combined file
    for seq in sequences:
        individual_gff_file = os.path.join(args.output, f"{seq.id}_fimo.gff")
        if os.path.exists(individual_gff_file):
            os.remove(individual_gff_file)

    # Remove fimo_output and graphics folders from the main output directory
    fimo_output_path = os.path.join(args.output, "fimo_out")
    output_graphics = os.path.join(args.output, "graphics")
    if os.path.exists(fimo_output_path):
        shutil.rmtree(fimo_output_path)
    if os.path.exists(output_graphics):
        shutil.rmtree(output_graphics)        

    # Write the promoter sequences to the output FASTA file
    output_fasta = os.path.join(args.output, "motif_sequences.fasta")
    os.makedirs(args.output, exist_ok=True)
    write_promoter_fasta(promoter_sequences, output_fasta)

    # Draw sequence graphics
    output_graphics = os.path.join(args.output, "output_graphics.png")
    output_pdf = os.path.join(args.output, "output_graphics.pdf")
    no_motif_file = os.path.join(args.output, "seqs_no_motifs.list")
    draw_sequence_graphics(sequences, all_promoter_regions, output_graphics, output_pdf, no_motif_file)
