import os
import re
import sys
import argparse
import subprocess
import tempfile
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_fasta(file_path):
    with open(file_path, 'r') as fasta_file:
        return list(SeqIO.parse(fasta_file, 'fasta'))

def scan_domains(sequence, hmm_file, threads, evalue):
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_fasta:
        SeqIO.write(sequence, temp_fasta, "fasta")
        temp_fasta.flush()

        output_file = tempfile.NamedTemporaryFile(mode="r", delete=False)
        log_file = os.path.join(args.output, "hmmscan.log")

        with open(log_file, "a") as log:  # Changed "w" to "a" to append to the log file
            print(f"Running hmmscan on sequence {sequence.id}")
            subprocess.run([
                "hmmscan",
                "--tblout", output_file.name,
                "-E", str(evalue),
                "--cpu", str(threads),
                hmm_file,
                temp_fasta.name
            ], stdout=log, stderr=subprocess.STDOUT)

        domain_regions = []

        with open(output_file.name, "r") as tblout:
            for line in tblout:
                if not line.startswith("#"):
                    match = re.search(r'(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)', line)
                    if match:
                        domain_name = match.group(1)
                        start, end = int(match.group(2)), int(match.group(3))
                        domain_regions.append((start, end, domain_name))

        os.remove(temp_fasta.name)
        os.remove(output_file.name)

    return domain_regions



def draw_sequence_graphics(sequences, domain_regions_list, output_path, output_pdf):
    fig, axes = plt.subplots(len(sequences), 1, figsize=(10, len(sequences) * 2))
    if len(sequences) == 1:
        axes = [axes]

    for ax, sequence, domain_regions in zip(axes, sequences, domain_regions_list):
        if domain_regions:
            for i, (start, end, domain_name) in enumerate(domain_regions):
                width = end - start + 1
                label_text = f"{domain_name}: {start}-{end}"
                rect = patches.Rectangle((start, 0), width, 1, linewidth=1, edgecolor='none', facecolor=plt.cm.viridis(float(i) / len(domain_regions)), alpha=0.5, label=label_text)
                ax.add_patch(rect)
            ax.legend(loc='upper right', fontsize='small')
        else:
            ax.text(0.5, 0.5, "No domains found", fontsize=12, ha='center', va='center', transform=ax.transAxes)

        ax.set_xlim(0, len(sequence))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_title(sequence.id)

    plt.tight_layout()
    plt.savefig(output_path)

    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig)

    plt.close()



def write_output_fasta(sequences, domain_regions_list, output_path):
    output_sequences = []

    for sequence, domain_regions in zip(sequences, domain_regions_list):
        domain_names = [domain_name for _, _, domain_name in domain_regions]
        domain_string = ";".join(domain_names)
        new_id = f"{sequence.id}_domains__{domain_string}"
        output_seq = SeqRecord(sequence.seq, id=new_id, description="")
        output_sequences.append(output_seq)

    with open(output_path, "w") as output_file:
        SeqIO.write(output_sequences, output_file, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output directory", required=True)
    parser.add_argument("-hmm", "--hmm", help="Path to the HMM file", required=True)
    parser.add_argument("-e", "--evalue", help="E-value threshold for hmmscan", default=10.0, type=float)
    parser.add_argument("-t", "--threads", help="Number of threads to use for hmmscan", default=1, type=int)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    sequences = parse_fasta(args.input)

    all_domain_regions = []
    for seq in sequences:
        domain_regions = scan_domains(seq, args.hmm, args.threads, args.evalue)
        all_domain_regions.append(domain_regions)

    output_graphics = os.path.join(args.output, "output_graphics.png")
    output_pdf = os.path.join(args.output, "output_graphics.pdf")
    draw_sequence_graphics(sequences, all_domain_regions, output_graphics, output_pdf)

    output_fasta = os.path.join(args.output, "output_sequences.fasta")
    write_output_fasta(sequences, all_domain_regions, output_fasta)
