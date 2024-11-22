import argparse
import csv
import os
from Bio import Entrez

# Configure email for Entrez
Entrez.email = "your_email@example.com"

def validate_accession(accession, db_type="nucleotide"):
    """Check if an accession is valid and get its summary."""
    try:
        handle = Entrez.esummary(db=db_type, id=accession, retmode="xml")
        summary = Entrez.read(handle)
        handle.close()
        return summary
    except Exception as e:
        print(f"Error validating accession {accession}: {e}")
        return None

def fetch_suppressed_replacement(accession, db_type="protein"):
    """Fetch a replacement accession for a suppressed record."""
    summary = validate_accession(accession, db_type)
    if not summary:
        return None

    for item in summary:
        if "AccessionVersion" in item:
            return item["AccessionVersion"]

    print(f"No replacement accession found for {accession}.")
    return None

def resolve_to_accession(gi_id, db_type="nucleotide"):
    """Resolve a GI number to a proper accession number."""
    try:
        summary = validate_accession(gi_id, db_type)
        if summary:
            for item in summary:
                if "AccessionVersion" in item:
                    return item["AccessionVersion"]
    except Exception as e:
        print(f"Error resolving GI number {gi_id} to accession: {e}")
    return gi_id

def fetch_cross_references(accession, from_db, to_db):
    """Find linked accession using Entrez elink API and resolve GI numbers to accessions."""
    try:
        handle = Entrez.elink(dbfrom=from_db, db=to_db, id=accession)
        records = Entrez.read(handle)
        handle.close()

        cross_refs = []
        for record in records:
            if "LinkSetDb" in record:
                for linkset in record["LinkSetDb"]:
                    for link in linkset["Link"]:
                        linked_id = link["Id"]
                        resolved_id = resolve_to_accession(linked_id, db_type=to_db)
                        cross_refs.append(resolved_id)
        return cross_refs
    except Exception as e:
        print(f"Error linking '{accession}' from '{from_db}' to '{to_db}': {e}")
        return []

def fetch_fasta(accession, db_type="nucleotide"):
    """Fetch the FASTA sequence for an accession."""
    try:
        handle = Entrez.efetch(db=db_type, id=accession, rettype="fasta", retmode="text")
        fasta = handle.read()
        handle.close()

        lines = fasta.splitlines()
        header = lines[0]
        sequence = "".join(lines[1:])  # Join sequence lines into a single string
        if len(sequence) > 10000:
            print(f"Sequence for {accession} is too long. Omitting.")
            return None

        return f"{header}\n{sequence}"
    except Exception as e:
        print(f"Error fetching FASTA for accession {accession}: {e}")
        return None

def construct_fasta_header(species, gene_name, alt_gene_name):
    """Construct the FASTA header as per the specifications."""
    species_parts = species.split("_")
    species_prefix = species_parts[0][:3] + species_parts[1][:3]
    return f">{species_prefix}_{gene_name}_{alt_gene_name}"

def fetch_and_add_missing_data(row, field, species, gene_name, alt_gene_name, fasta_files):
    """Fetch missing accession numbers and their FASTAs."""
    accession = row.get(f"{field}_Accession")

    if accession:
        fasta = fetch_fasta(accession, db_type="nucleotide" if field != "Protein" else "protein")
        if not fasta and field == "Protein":
            replacement = fetch_suppressed_replacement(accession, db_type="protein")
            if replacement:
                print(f"Using replacement accession {replacement} for {accession}.")
                fasta = fetch_fasta(replacement, db_type="protein")
                row[f"{field}_Accession"] = replacement

        if fasta:
            header = construct_fasta_header(species, gene_name, alt_gene_name)
            sequence = "".join(fasta.splitlines()[1:])  # Single-line sequence
            row[f"{field}_FASTA"] = f"{header} {sequence}"
            fasta_files[field].write(f"{header}\n{sequence}\n")
        else:
            row[f"{field}_FASTA"] = f"NOTE: Sequence omitted for accession {accession} (too long or no region found)"
            print(f"Sequence omitted for {field} accession {accession}")
    else:
        print(f"No accession provided for {field}.")

def resolve_missing_ids(row):
    """Find and populate missing IDs using cross-references."""
    db_mappings = {
        "Protein_Accession": ("protein", [("nucleotide", "mRNA_Accession")]),
        "mRNA_Accession": ("nucleotide", [("protein", "Protein_Accession")]),
    }

    for current_field, (db_from, targets) in db_mappings.items():
        accession = row.get(current_field)
        if accession:
            for db_to, target_field in targets:
                if not row.get(target_field):
                    cross_refs = fetch_cross_references(accession, db_from, db_to)
                    if cross_refs:
                        row[target_field] = cross_refs[0]
                        print(f"Resolved {target_field} from {current_field} using accession {accession}.")
                    else:
                        print(f"No cross-references found for {current_field} ({accession}) in {db_to}.")

def process_tsv(input_file, output_file):
    """Process the TSV file."""
    basename = os.path.splitext(output_file)[0]
    mrna_fasta_file = f"{basename}.mrna.fa"
    protein_fasta_file = f"{basename}.prot.fa"

    with open(mrna_fasta_file, "w") as mrna_file, open(protein_fasta_file, "w") as protein_file:
        fasta_files = {"mRNA": mrna_file, "Protein": protein_file}

        with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
            reader = csv.DictReader(infile, delimiter="\t")
            fieldnames = [field for field in reader.fieldnames if not field.startswith("Gene_")]
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()

            for row in reader:
                species = row["Species"]
                gene_name = row["Gene_name"]
                alt_gene_name = row["Alternate_Gene_Name"]

                resolve_missing_ids(row)

                for field in ["mRNA", "Protein"]:
                    if not row.get(f"{field}_FASTA"):
                        fetch_and_add_missing_data(row, field, species, gene_name, alt_gene_name, fasta_files)

                writer.writerow({key: row[key] for key in fieldnames})

def main():
    parser = argparse.ArgumentParser(description="Complete a TSV file with missing accession & FASTA data.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    process_tsv(args.input, args.output)

if __name__ == "__main__":
    main()
