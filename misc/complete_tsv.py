import argparse
import csv
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

def fetch_cross_references(accession, from_db, to_db):
    """Find linked accession using Entrez elink API."""
    try:
        handle = Entrez.elink(dbfrom=from_db, db=to_db, id=accession)
        records = Entrez.read(handle)
        handle.close()

        # Extract linked accessions
        cross_refs = []
        for record in records:
            if "LinkSetDb" in record:
                for linkset in record["LinkSetDb"]:
                    for link in linkset["Link"]:
                        cross_refs.append(link["Id"])
        return cross_refs
    except Exception as e:
        print(f"Error linking '{accession}' from '{from_db}' to '{to_db}': {e}")
        return []

def resolve_missing_ids(row):
    """Find and populate missing IDs using cross-references."""
    db_mappings = {
        "Protein_Accession": ("protein", [("nucleotide", "Gene_Accession"), ("nucleotide", "mRNA_Accession")]),
        "Gene_Accession": ("nucleotide", [("protein", "Protein_Accession"), ("nucleotide", "mRNA_Accession")]),
        "mRNA_Accession": ("nucleotide", [("protein", "Protein_Accession"), ("nucleotide", "Gene_Accession")]),
    }

    for current_field, (db_from, targets) in db_mappings.items():
        accession = row.get(current_field)
        if accession:  # Proceed only if the current field has an accession
            for db_to, target_field in targets:
                if not row.get(target_field):  # Populate only if target field is missing
                    cross_refs = fetch_cross_references(accession, db_from, db_to)
                    if cross_refs:
                        row[target_field] = cross_refs[0]  # Use the first linked accession
                        print(f"Resolved {target_field} from {current_field} using accession {accession}.")
                    else:
                        print(f"No cross-references found for {current_field} ({accession}) in {db_to}.")

def fetch_fasta(accession, db_type="nucleotide"):
    """Fetch the FASTA sequence for an accession."""
    try:
        handle = Entrez.efetch(db=db_type, id=accession, rettype="fasta", retmode="text")
        fasta = handle.read()
        handle.close()
        return fasta
    except Exception as e:
        print(f"Error fetching FASTA for accession {accession}: {e}")
        return None

def construct_fasta_header(species, gene_name, alt_gene_name):
    """Construct the FASTA header as per the specifications."""
    species_parts = species.split("_")
    species_prefix = species_parts[0][:3] + species_parts[1][:3]
    return f">{species_prefix}_{gene_name}_{alt_gene_name}\n"

def fetch_and_add_missing_data(row, field, species, gene_name, alt_gene_name):
    """Fetch missing accession numbers and their FASTAs."""
    accession = row.get(f"{field}_Accession")

    if accession:
        fasta = fetch_fasta(accession, db_type="nucleotide" if field != "Protein" else "protein")
        if not fasta and field == "Protein":  # Attempt to resolve suppressed accessions for proteins
            replacement = fetch_suppressed_replacement(accession, db_type="protein")
            if replacement:
                print(f"Using replacement accession {replacement} for {accession}.")
                fasta = fetch_fasta(replacement, db_type="protein")
                row[f"{field}_Accession"] = replacement

        if fasta:
            header = construct_fasta_header(species, gene_name, alt_gene_name)
            row[f"{field}_FASTA"] = header + "".join(fasta.split("\n")[1:])
        else:
            print(f"Error fetching FASTA for {field} accession {accession}")
    else:
        print(f"No accession provided for {field}.")

def process_tsv(input_file, output_file):
    """Process the TSV file."""
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        
        # Ensure FASTA fields are added only once
        fasta_fields = ["Gene_FASTA", "mRNA_FASTA", "Protein_FASTA"]
        fieldnames = reader.fieldnames + [field for field in fasta_fields if field not in reader.fieldnames]
        
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            species = row["Species"]
            gene_name = row["Gene_name"]
            alt_gene_name = row["Alternate_Gene_Name"]

            # Attempt to resolve missing IDs
            resolve_missing_ids(row)

            # Fetch FASTA sequences for each field
            for field in ["Gene", "mRNA", "Protein"]:
                if not row.get(f"{field}_FASTA"):
                    fetch_and_add_missing_data(row, field, species, gene_name, alt_gene_name)

            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Complete a TSV file with missing accession FASTA data.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    process_tsv(args.input, args.output)

if __name__ == "__main__":
    main()
