import argparse
import logging
import os
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pipeline_utils

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser(description="Pipeline for extracting gene regions, upstream sequences, and promoter prediction.")
    parser.add_argument("-m", "--mrna_fasta", required=True, help="mRNA sequences in FASTA format.")
    parser.add_argument("-g", "--genome_fasta", required=True, help="Genome sequence in FASTA format.")
    parser.add_argument("-u", "--upstream_length", type=int, required=True, help="Length of nucleotides upstream of the start codon.")
    parser.add_argument("-o", "--output_dir", default=".", help="Output directory for the resulting FASTA file.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for parallel processing.")
    parser.add_argument("-e", "--evalue", type=float, default=0.001, help="E-value threshold for BLAST search.")
    parser.add_argument("-p", "--promoter", action='store_true', help="Run promoter prediction on the extracted sequences.")
    args = parser.parse_args()

    # Create output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse input files
    try:
        mrna_sequences, genome_sequences = pipeline_utils.parse_input_files(args.mrna_fasta, args.genome_fasta)
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        return

    # Create BLAST database using subprocess
    genome_db = pipeline_utils.create_blast_database(args.genome_fasta)

    # Perform BLAST search using subprocess
    blast_output_file = pipeline_utils.perform_blast_search(args.mrna_fasta, genome_db, args.evalue, args.threads, args.output_dir)

    # Parse BLAST output
    with open(blast_output_file) as blast_output_handle:
        blast_records = list(NCBIXML.parse(blast_output_handle))

    # Extract upstream sequences
    extracted_sequences = pipeline_utils.extract_upstream_sequences(blast_records, genome_sequences, args.upstream_length)
    output_fasta_file = os.path.join(args.output_dir, "extracted_sequences.fasta")
    pipeline_utils.output_fasta(extracted_sequences, output_fasta_file)

    # Perform promoter prediction if required
    if args.promoter:
        promoter_output_file = os.path.join(args.output_dir, "promoter_predictions.txt")
        pipeline_utils.predict_promoters(output_fasta_file, promoter_output_file)

if __name__ == "__main__":
    main()
