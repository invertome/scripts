import argparse
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio.Blast import NCBIXML

def parse_input_files(mrna_fasta, genome_fasta):
    """
    Read and parse input files using Biopython's SeqIO.parse.
    """
    mrna_sequences = list(SeqIO.parse(mrna_fasta, "fasta"))
    genome_sequence = next(SeqIO.parse(genome_fasta, "fasta"))

    return mrna_sequences, genome_sequence

def perform_blast_search(mrna_sequences, genome_fasta, blast_output_file):
    """
    Perform BLAST search using NcbiblastnCommandline.
    """
    # Create BLAST database for the genome sequence with parsed sequence IDs
    makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=genome_fasta, parse_seqids=True)
    makeblastdb_cline()

    # Perform BLAST search for each mRNA sequence
    blast_outputs = []
    for mrna in mrna_sequences:
        # Create a temporary FASTA file for the mRNA sequence
        mrna_fasta = f"temp_mrna_{mrna.id}.fasta"
        SeqIO.write(mrna, mrna_fasta, "fasta")

        # Run BLAST
        blastn_cline = NcbiblastnCommandline(query=mrna_fasta, db=genome_fasta, outfmt=5, out=blast_output_file)
        blastn_cline()

        # Parse BLAST results
        blast_records = list(NCBIXML.parse(open(blast_output_file)))
        blast_outputs.append(blast_records)

        # Remove temporary files
        os.remove(mrna_fasta)

    return blast_outputs

def extract_upstream_sequences(blast_outputs, genome_sequence, upstream_length):
    """
    Extract the specified length of nucleotides upstream of the start codon.
    """
    extracted_sequences = []

    for blast_records in blast_outputs:
        for alignment in blast_records[0].alignments:
            for hsp in alignment.hsps:
                start = max(0, hsp.sbjct_start - upstream_length - 1)
                end = hsp.sbjct_start - 1
                extracted_seq = genome_sequence.seq[start:end]
                extracted_sequences.append(SeqIO.SeqRecord(extracted_seq, id=alignment.title))

    return extracted_sequences

def output_fasta(sequences, output_filename):
    """
    Output extracted sequences as FASTA using SeqIO.write.
    """
    SeqIO.write(sequences, output_filename, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Pipeline for extracting gene regions and upstream sequences.")
    parser.add_argument("-m", "--mrna_fasta", required=True, help="mRNA sequences in FASTA format.")
    parser.add_argument("-g", "--genome_fasta", required=True, help="Genome sequence in FASTA format.")
    parser.add_argument("-u", "--upstream_length", type=int, required=True, help="Length of nucleotides upstream of the start codon.")
    parser.add_argument("-o", "--output_dir", default=".", help="Output directory for the resulting FASTA file.")
    args = parser.parse_args()

    # Create output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Run the pipeline
    mrna_sequences, genome_sequence = parse_input_files(args.mrna_fasta, args.genome_fasta)
    blast_output_file = os.path.join(args.output_dir, "blast_output.xml")
    blast_outputs = perform_blast_search(mrna_sequences, args.genome_fasta, blast_output_file)
    extracted_sequences = extract_upstream_sequences(blast_outputs, genome_sequence, args.upstream_length)
    output_fasta_file = os.path.join(args.output_dir, "extracted_sequences.fasta")
    output_fasta(extracted_sequences, output_fasta_file)

if __name__ == "__main__":
    main()
