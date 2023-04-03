import os
import sys
import argparse
import subprocess
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_fasta(file_path):
    with open(file_path, 'r') as fasta_file:
        return list(SeqIO.parse(fasta_file, 'fasta'))

def scan_domains(sequence, hmm_file):
    with open("temp_sequence.fasta", "w") as temp_file:
        SeqIO.write(sequence, temp_file, "fasta")

    output_path = os.path.join(args.output, "hmmer_out")
    os.makedirs(output_path, exist_ok=True)

    subprocess.run([
        "hmmsearch",
        "--tblout", os.path.join(output_path, f"{sequence.id}.tblout"),
        hmm_file,
        "temp_sequence.fasta"
    ])

    os.remove("temp_sequence.fasta")

    domains = []

    with open(os.path.join(output_path, f"{sequence.id}.tblout"), "r") as tblout_file:
        for line in tblout_file:
            if not line.startswith("#"):
                columns = line.strip().split()
                if len(columns) >= 19:
                    target_name, accession, tlen, query_name, accession2, qlen, seq_evalue, seq_score, seq_bias, domain_num, domain_total, c_evalue, i_evalue, dom_score, dom_bias, hmm_start, hmm_end, ali_start, ali_end = columns[:19]
                    domains.append((int(ali_start) - 1, int(ali_end), query_name))

    return domains

def draw_sequence_graphics(sequences, domain_regions_list, output_path, output_pdf):
    fig, axes = plt.subplots(len(sequences), 1, figsize=(10, len(sequences) * 2))
    if len(sequences) == 1:
        axes = [axes]

    for ax, sequence, domain_regions in zip(axes, sequences, domain_regions_list):
        if domain_regions:
            for i, (start, end, query_name) in enumerate(domain_regions):
                label_text = f"Domain {i+1}: {start}-{end}, {query_name}"
                ax.axvspan(start, end, color=plt.cm.viridis(float(i) / len(domain_regions)), alpha=0.5, label=label_text)
            ax.legend(loc='upper right', fontsize='small')

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
    annotated_sequences = []

    for sequence, domain_regions in zip(sequences, domain_regions_list):
        if domain_regions:
            domain_names = [query_name for _, _, query_name in domain_regions]
            sequence.id = f"{sequence.id} domains: {', '.join(domain_names)}"
        annotated_sequences.append(sequence)

    with open(output_path, "w") as output_file:
        for seq in annotated_sequences:
            SeqIO.write(seq, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to the input FASTA file", required=True)
    parser.add_argument("-o", "--output", help="Path to the output directory", required=True)
    parser.add_argument("-m", "--hmm", help="Path to the HMM file", required=True)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    sequences = parse_fasta(args.input)

    all_domain_regions = []
    for seq in sequences:
        domain_regions = scan_domains(seq, args.hmm)
        all_domain_regions.append(domain_regions)

    output_graphics = os.path.join(args.output, "output_graphics.png")
    output_pdf = os.path.join(args.output, "output_graphics.pdf")
    draw_sequence_graphics(sequences, all_domain_regions, output_graphics, output_pdf)

    output_fasta = os.path.join(args.output, "output_sequences.fasta")
    write_output_fasta(sequences, all_domain_regions, output_fasta)

