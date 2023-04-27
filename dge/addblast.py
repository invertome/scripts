import csv
import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Add protein IDs and percent identity to gene table CSV')
parser.add_argument('-i', '--input', nargs='+', required=True, help='Input CSV files (space-separated)')
parser.add_argument('-b', '--blast', required=True, help='BLAST output file')
args = parser.parse_args()

# Read the BLAST output and store the results in a dictionary
blast_dict = {}
with open(args.blast, 'r') as blast_file:
    for line in blast_file:
        fields = line.strip().split('\t')
        gene_id = fields[0]
        protein_id = fields[1]
        percent_identity = fields[2]
        blast_dict[gene_id] = {'protein_id': protein_id, 'percent_identity': percent_identity}

# Process each input file
for input_file in args.input:
    # Generate output filename from input file basename
    input_basename, _ = os.path.splitext(input_file)
    output_file = f'{input_basename}_annotated.csv'

    # Read the CSV file and add new columns for protein IDs and percent identity
    with open(input_file, 'r') as csvfile_in, open(output_file, 'w', newline='') as csvfile_out:
        csv_reader = csv.reader(csvfile_in)
        csv_writer = csv.writer(csvfile_out)

        # Write the header row
        header = next(csv_reader)
        header.extend(['ProteinID', 'PercentIdentity'])
        csv_writer.writerow(header)

        # Write the data rows with the new columns for protein IDs and percent identity
        for row in csv_reader:
            gene_id = row[0]
            protein_data = blast_dict.get(gene_id, {'protein_id': 'NA', 'percent_identity': 'NA'})
            row.extend([protein_data['protein_id'], protein_data['percent_identity']])
            csv_writer.writerow(row)
