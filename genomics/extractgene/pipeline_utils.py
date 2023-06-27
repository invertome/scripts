import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio.Blast import NCBIXML
import subprocess

def parse_input_files(mrna_fasta, genome_fasta):
    """
    Read and parse input files using Biopython's SeqIO.parse.
    """
    mrna_sequences = list(SeqIO.parse(mrna_fasta, "fasta"))
    genome_sequence = next(SeqIO.parse(genome_fasta, "fasta"))

    return mrna_sequences, genome_sequence

def perform_blast_search(mrna, genome_fasta, blast_output_file, evalue):
    """
    Perform BLAST search using NcbiblastnCommandline.
    """
    # Create BLAST database for the genome sequence with parsed sequence IDs
    makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=genome_fasta, parse_seqids=True)
    makeblastdb_cline()

    # Create a temporary FASTA file for the current mRNA sequence
    temp_mrna_fasta = os.path.join(os.path.dirname(genome_fasta), f"{mrna.id}_temp.fasta")
    SeqIO.write(mrna, temp_mrna_fasta, "fasta")

    # Perform the BLAST search
    blastn_cline = NcbiblastnCommandline(query=temp_mrna_fasta, db=genome_fasta, evalue=evalue, outfmt=5, out=blast_output_file)
    blastn_cline()

    # Remove the temporary FASTA file
    os.remove(temp_mrna_fasta)

    # Parse the BLAST output
    with open(blast_output_file, "r") as blast_output_handle:
        blast_record = NCBIXML.read(blast_output_handle)

    return blast_record

def extract_upstream_sequences(blast_outputs, genome_sequence, upstream_length):
    extracted_sequences = []

    for i, blast_record in enumerate(blast_outputs):
        if blast_record.alignments:
            alignment = blast_record.alignments[0]  # only take top alignment
            if alignment.hsps:
                hsp = alignment.hsps[0]  # only take top hsp

                start = hsp.sbjct_start - upstream_length
                end = hsp.sbjct_start - 1

                if start < 1:
                    start = 1

                if start < end and end <= len(genome_sequence):  # check coordinates are valid
                    extracted_seq = genome_sequence.seq[start - 1:end]
                    print(f"Extracted sequence for {blast_record.query}: {extracted_seq}")  # Debug print statement
                    extracted_sequences.append((blast_record.query, extracted_seq))  # Keep original ID
        else:
            print(f"No alignments found for {blast_record.query}")  # Debug print statement

    return extracted_sequences



def output_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        for i, (original_id, seq) in enumerate(sequences):
            output_handle.write(f">{original_id}_seq_{i}\n{seq}\n")  # Include original ID in header

def predict_promoters(input_fasta_file, output_file):
    """
    Run promoter prediction using the Promoter2.0 tool.
    """
    promoter_cline = f"promoter2 {input_fasta_file} > {output_file}"
    subprocess.run(promoter_cline, shell=True, check=True)
