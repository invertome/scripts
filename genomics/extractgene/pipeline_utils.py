import os
from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
import logging
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages


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

def create_hisat2_index(genome_fasta):
    index_base = os.path.splitext(genome_fasta)[0]
    index_files = [f"{index_base}.{i}.ht2" for i in range(1, 9)]

    # Check if index files already exist
    if all(os.path.exists(f) for f in index_files):
        logging.info(f"HISAT2 index files already exist for {genome_fasta}. Skipping index creation.")
    else:
        logging.info(f"Creating HISAT2 index for {genome_fasta}.")
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

def extract_entire_gene(sam_file, genome_sequences):
    gene_sequences = {}
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.split('\t')
            query_name = fields[0]
            ref_name = fields[2]
            pos = int(fields[3])
            cigar = fields[5]

            if ref_name in genome_sequences:
                ref_seq = genome_sequences[ref_name]
                if query_name not in gene_sequences:
                    gene_sequences[query_name] = {"start": pos, "end": pos, "exons": []}
                
                current_pos = pos
                # Track aligned regions based on CIGAR string
                for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
                    length = int(length)
                    if op in "M=X":  # Match or alignment match
                        gene_sequences[query_name]["exons"].append((current_pos, current_pos + length))
                        current_pos += length
                    elif op == 'D':  # Deletion in the reference
                        current_pos += length
                    elif op == 'N':  # Skipped region (intron)
                        current_pos += length
                    elif op in 'IS':  # Insertion to the reference or soft clipping, skip in reference
                        continue
                
                gene_sequences[query_name]["start"] = min(gene_sequences[query_name]["start"], pos)
                gene_sequences[query_name]["end"] = max(gene_sequences[query_name]["end"], current_pos)

    extracted_genes = []
    for gene, coords in gene_sequences.items():
        start = coords["start"]
        end = coords["end"]
        seq = ref_seq.seq[start-1:end]

        # Constructing the full gene sequence with aligned regions in uppercase and non-aligned (introns) in lowercase
        gene_seq = ''
        last_pos = start - 1
        for exon_start, exon_end in sorted(coords["exons"]):
            # Add non-aligned (intron) sequence
            if exon_start > last_pos:
                gene_seq += str(seq[last_pos-start:exon_start-start].lower())
            # Add aligned (exon) sequence
            gene_seq += str(seq[exon_start-start:exon_end-start].upper())
            last_pos = exon_end
        
        # Add any trailing non-aligned sequence
        if last_pos < end:
            gene_seq += str(seq[last_pos-start:].lower())
        extracted_genes.append((gene, gene_seq, coords["exons"]))

    return extracted_genes

def output_gene_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        for gene, seq, _ in sequences:
            output_handle.write(f">{gene}_full_gene\n{seq}\n")

def draw_intron_exon_structure(gene_sequences, output_pdf):
    with PdfPages(output_pdf) as pdf:
        for gene, _, exons in gene_sequences:
            fig, ax = plt.subplots(figsize=(10, 2))

            gene_start = min([start for start, _ in exons])
            gene_end = max([end for _, end in exons])
            ax.plot([gene_start, gene_end], [1, 1], color='black', lw=1)

            for start, end in exons:
                if exons.index((start, end)) > 0:
                    previous_end = exons[exons.index((start, end)) - 1][1]
                    if start > previous_end:
                        ax.plot([previous_end, start], [1, 1], color='red', lw=1)
                ax.add_patch(patches.Rectangle((start, 0.9), end-start, 0.2, edgecolor='black', facecolor='blue'))

            ax.text(gene_start, 0.8, f"{gene_start}", fontsize=8, ha='center')
            ax.text(gene_end, 0.8, f"{gene_end}", fontsize=8, ha='center')

            ax.set_ylim(0.7, 1.3)
            ax.set_xlim(gene_start - 100, gene_end + 100)
            ax.axis('off')
            ax.set_title(f"{gene} (Genomic Positions: {gene_start}-{gene_end})")

            pdf.savefig(fig)
            plt.close(fig)
