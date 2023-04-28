import argparse
import subprocess
import os
from Bio import SeqIO

# Replaces IDs in the results file with the species names
def replace_ids(species_tab_file, results_file, output_file):
    id_to_species = {}
    with open(species_tab_file, 'r') as f:
        for line in f:
            _, id_, species, *rest = line.split()
            id_to_species[id_] = species

    with open(results_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                columns = line.split()
                target_id = columns[-1]
                species = id_to_species.get(target_id, target_id)
                new_line = "\t".join(columns[:-1] + [species]) + "\n"
                outfile.write(new_line)
            else:
                outfile.write(line)

# Extracts sequences from the database (or FASTA file) and writes them to an output FASTA file with updated headers
def extract_fasta_sequences(database, fasta_file, species_tab_file, results_file, output_fasta_file):
    # Collect target IDs from the results file
    target_ids = set()
    with open(results_file, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                target_id = line.split()[0]
                target_ids.add(target_id)

    # Create a dictionary mapping sequence IDs to species names
    id_to_species = {}
    with open(species_tab_file, 'r') as f:
        for line in f:
            _, id_, species, *rest = line.split()
            id_to_species[id_] = species

    if not database:
        # Create BLAST database if not provided
        database = "temp_blast_db"
        subprocess.run(["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", database])

    with open(output_fasta_file, 'w') as out_fasta:
        for target_id in target_ids:
            # Retrieve sequence from BLAST database
            cmd = ["blastdbcmd", "-db", database, "-entry", target_id, "-outfmt", "%f"]
            sequence = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout

            # Parse the retrieved sequence
            record = SeqIO.read(io.StringIO(sequence), "fasta")

            # Update the header with the species ID
            species = id_to_species.get(record.id, record.id)
            record.description = f"{record.description} {species}"

            # Write the updated sequence to the output FASTA file
            SeqIO.write(record, out_fasta, "fasta")

    if not database:
        # Remove temporary BLAST database files
        for ext in [".nhr", ".nin", ".nsq"]:
            try:
                os.remove(f"{database}{ext}")
            except FileNotFoundError:
                pass

def main():
    parser = argparse.ArgumentParser(description='Process OrthoDB HMM output')
    parser.add_argument('-s', '--species-tab', required=True, help='Path to the species tab file')
    parser.add_argument('-r', '--results', required=True, help='Path to the results file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')
    parser.add_argument('-f', '--fasta', help='Path to the input FASTA file (optional)')
    parser.add_argument('-d', '--db', help='Path to the existing BLAST database (optional)')
    args = parser.parse_args()

    # Replace IDs with species names in the results file
    replace_ids(args.species_tab, args.results, args.output)

    if args.fasta or args.db:
        # Extract sequences from the database (or FASTA file) and write them to an output FASTA file with updated headers
        output_fasta_file = args.output.rsplit('.', 1)[0] + '_sequences.fasta'
        extract_fasta_sequences(args.db, args.fasta, args.species_tab, args.results, output_fasta_file)

if __name__ == '__main__':
    main()
