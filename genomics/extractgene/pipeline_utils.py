import os
from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_input_files(mrna_fasta, genome_fasta):
    mrna_sequences = list(SeqIO.parse(mrna_fasta, "fasta"))
    genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    return mrna_sequences, genome_sequences

def create_blast_database(genome_fasta):
    db_name = os.path.splitext(genome_fasta)[0]
    makeblastdb_cline = ['makeblastdb', '-in', genome_fasta, '-dbtype', 'nucl', '-out', db_name]
    subprocess.run(makeblastdb_cline, check=True)
    return db_name

def perform_blast_search(mrna_fasta, genome_db, evalue, threads, output_dir):
    blast_output_file = os.path.join(output_dir, "blast_results.xml")
    blastn_cline = [
        'blastn', '-query', mrna_fasta, '-db', genome_db, '-out', blast_output_file,
        '-outfmt', '5', '-evalue', str(evalue), '-num_threads', str(threads)
    ]
    subprocess.run(blastn_cline, check=True)
    return blast_output_file

def extract_upstream_sequences(blast_records, genome_sequences, upstream_length):
    extracted_sequences = []
    for blast_record in blast_records:
        if blast_record.alignments:
            alignment = blast_record.alignments[0]  # only take top alignment
            if alignment.hsps:
                hsp = alignment.hsps[0]  # only take top hsp
                reference_id = alignment.hit_def.split()[0]  # split at the first whitespace

                if reference_id not in genome_sequences:
                    logging.error(f"Reference {reference_id} not found in genome sequences.")
                    continue

                ref_seq = genome_sequences[reference_id]
                start = hsp.sbjct_start - upstream_length
                end = hsp.sbjct_start - 1

                if start < 1:
                    start = 1

                if start < end and end <= len(ref_seq):  # check coordinates are valid
                    extracted_seq = ref_seq.seq[start - 1:end]
                    logging.info(f"Extracted sequence for {blast_record.query}: {extracted_seq}")
                    extracted_sequences.append((blast_record.query, extracted_seq))  # Keep original ID
                else:
                    logging.error(f"Invalid coordinates for {blast_record.query}: start={start}, end={end}, len={len(ref_seq)}")
        else:
            logging.warning(f"No alignments found for {blast_record.query}")
    return extracted_sequences

def output_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        for i, (original_id, seq) in enumerate(sequences):
            output_handle.write(f">{original_id}_seq_{i}\n{seq}\n")

def predict_promoters(input_fasta_file, output_file):
    promoter_cline = ["promoter2", input_fasta_file, output_file]
    subprocess.run(promoter_cline, check=True)
