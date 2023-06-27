import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio.Blast import NCBIXML
import subprocess
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_input_files(mrna_fasta, genome_fasta):
    """
    Read and parse input files using Biopython's SeqIO.parse.
    """
    mrna_sequences = list(SeqIO.parse(mrna_fasta, "fasta"))
    genome_sequence = next(SeqIO.parse(genome_fasta, "fasta"))

    return mrna_sequences, genome_sequence

def perform_blast_search(mrna, genome_db, evalue):
    """
    Perform BLAST search using NcbiblastnCommandline.
    """
    # Create a temporary FASTA file for the current mRNA sequence
    temp_mrna_fasta = Path(genome_db.parent, f"{mrna.id}_temp.fasta")
    SeqIO.write(mrna, temp_mrna_fasta, "fasta")

    # Temporary output file
    blast_output_file = Path(genome_db.parent, f"{mrna.id}_blast.xml")
    
    # Perform the BLAST search
    blastn_cline = NcbiblastnCommandline(query=str(temp_mrna_fasta), db=str(genome_db), evalue=evalue, outfmt=5, out=str(blast_output_file))
    blastn_cline()

    # Parse the BLAST output
    with blast_output_file.open() as blast_output_handle:
        blast_record = NCBIXML.read(blast_output_handle)

    # Clean up temporary files
    temp_mrna_fasta.unlink()
    blast_output_file.unlink()

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
                    logging.info(f"Extracted sequence for {blast_record.query}: {extracted_seq}")  # Debug print statement
                    extracted_sequences.append((blast_record.query, extracted_seq))  # Keep original ID
                else:
                    logging.error(f"Invalid coordinates for {blast_record.query}: start={start}, end={end}, len={len(genome_sequence)}")
        else:
            logging.warning(f"No alignments found for {blast_record.query}")  # Debug print statement

    return extracted_sequences

def output_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        for i, (original_id, seq) in enumerate(sequences):
            output_handle.write(f">{original_id}_seq_{i}\n{seq}\n")  # Include original ID in header

def predict_promoters(input_fasta_file, output_file):
    """
    Run promoter prediction using the Promoter2.0 tool.
    """
    promoter_cline = ["promoter2", input_fasta_file, ">", output_file]
    subprocess.run(promoter_cline, shell=True, check=True)
