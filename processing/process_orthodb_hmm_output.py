import argparse
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO

# Check if the required package 'biopython' is installed
try:
    import Bio
except ImportError:
    print("Biopython package not found. Installing it now...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
    import Bio

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
def extract_fasta_sequences(fasta_file, species_tab_file, results_file, output_fasta_file, num_threads):
    target_ids = set()
    with open(results_file, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                target_id = line.split()[0]
                target_ids.add(target_id)

    id_to_species = {}
    with open(species_tab_file, 'r') as f:
        for line in f:
            _, id_, species, *rest = line.split()
            id_to_species[id_] = species.replace(" ", "_")

    indexed_fasta = SeqIO.index(fasta_file, "fasta")

    def write_sequence(target_id):
        record = indexed_fasta.get(target_id)
        if record is not None:
            species = id_to_species.get(target_id, target_id)
            record.description = f"{record.description} {species}"
            return record

    with open(output_fasta_file, 'w') as out_fasta:
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            records = list(executor.map(write_sequence, target_ids))

        for record in records:
            if record:
                SeqIO.write(record, out_fasta, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Process OrthoDB HMM output')
    parser.add_argument('-s', '--species-tab', required=True, help='Path to the species tab file')
    parser.add_argument('-r', '--results', required=True, help='Path to the results file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')
    parser.add_argument('-f', '--fasta', help='Path to the input FASTA file (optional)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads (default: 4)')
    args = parser.parse_args()

    if not args.species_tab or not args.results or not args.output:
        print("Error: Missing required arguments. Please check the usage and provide the necessary arguments.")
        sys.exit(1)

    if not args.fasta:
        print("Warning: No FASTA file provided. Only replacing IDs in the results file.")
    else:
        print("Processing FASTA file...")

    replace_ids(args.species_tab, args.results, args.output)

    if args.fasta:
        output_fasta_file = args.output.rsplit('.', 1)[0] + '_sequences.fasta'
        extract_fasta_sequences(args.fasta, args.species_tab, args.results, output_fasta_file, args.threads)

if __name__ == '__main__':
    main()
