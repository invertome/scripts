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

def extract_upstream_sequences_blast(blast_records, genome_sequences, upstream_length):
    extracted_sequences = []
    for blast_record in blast_records:
        if blast_record.alignments:
            first_hsp = None
            first_hsp_start = None
            alignment_id = None
            
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if first_hsp is None or hsp.sbjct_start < first_hsp_start:
                        first_hsp = hsp
                        first_hsp_start = hsp.sbjct_start
                        alignment_id = alignment.hit_def.split()[0]

            if first_hsp and alignment_id:
                if alignment_id not in genome_sequences:
                    logging.error(f"Reference {alignment_id} not found in genome sequences.")
                    continue

                ref_seq = genome_sequences[alignment_id]
                start = first_hsp.sbjct_start - upstream_length
                end = first_hsp.sbjct_start - 1

                if start < 1:
                    start = 1

                if start < end and end <= len(ref_seq):
                    extracted_seq = ref_seq.seq[start - 1:end]
                    logging.info(f"Extracted sequence for {blast_record.query}: {extracted_seq}")
                    extracted_sequences.append((blast_record.query, extracted_seq))
                else:
                    logging.error(f"Invalid coordinates for {blast_record.query}: start={start}, end={end}, len={len(ref_seq)}")
        else:
            logging.warning(f"No alignments found for {blast_record.query}")
    return extracted_sequences

def create_hisat2_index(genome_fasta):
    index_base = os.path.splitext(genome_fasta)[0]
    hisat2_build_cline = ['hisat2-build', genome_fasta, index_base]
    subprocess.run(hisat2_build_cline, check=True)
    return index_base

def perform_hisat2_mapping(mrna_fasta, genome_index, threads, output_dir, score_min="L,0,-4", pen_noncansplice=3):
    sam_output_file = os.path.join(output_dir, "hisat2_results.sam")
    hisat2_cline = [
        'hisat2', '-f', '-x', genome_index, '-U', mrna_fasta, '-S', sam_output_file,
        '--threads', str(threads),
        '--score-min', score_min,
        '--pen-noncansplice', str(pen_noncansplice)
    ]
    subprocess.run(hisat2_cline, check=True)
    return sam_output_file


def extract_upstream_sequences_hisat2(sam_file, genome_sequences, upstream_length):
    extracted_sequences = []
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.split('\t')
            query_name = fields[0]
            ref_name = fields[2]
            pos = int(fields[3])

            if ref_name in genome_sequences:
                ref_seq = genome_sequences[ref_name]
                start = pos - upstream_length
                end = pos - 1

                if start < 1:
                    start = 1

                if start < end and end <= len(ref_seq):
                    extracted_seq = ref_seq.seq[start - 1:end]
                    logging.info(f"Extracted sequence for {query_name}: {extracted_seq}")
                    extracted_sequences.append((query_name, extracted_seq))
                else:
                    logging.error(f"Invalid coordinates for {query_name}: start={start}, end={end}, len={len(ref_seq)}")
            else:
                logging.error(f"Reference {ref_name} not found in genome sequences.")

    return extracted_sequences

def output_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        for i, (original_id, seq) in enumerate(sequences):
            output_handle.write(f">{original_id}_seq_{i}\n{seq}\n")

def predict_promoters(input_fasta_file, output_file):
    promoter_cline = ["promoter2", input_fasta_file, output_file]
    subprocess.run(promoter_cline, check=True)
