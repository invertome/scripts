#!/bin/bash

set -e  # Exit on any error

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <taxon-name> <data-types> [--refseq]"
  echo "Example: $0 arthropoda genome protein rna --refseq"
  exit 1
fi

TAXON_NAME=$1
shift
REFSEQ_ONLY=false

# Check if the last argument is --refseq
if [ "${@: -1}" == "--refseq" ]; then
  REFSEQ_ONLY=true
  set -- "${@:1:$(($#-1))}"
fi

DATA_TYPES=$(echo $@ | tr ' ' ',')

# Set the number of cores to use
NUM_CORES=4  # Adjust this number based on your system's capability

# Function to download and hydrate RefSeq data
download_data() {
  if [ ! -d "${TAXON_NAME}_ncbi" ]; then
    echo "Downloading NCBI data for ${TAXON_NAME} with data types: ${DATA_TYPES}..."
    datasets download genome taxon "${TAXON_NAME}" --dehydrated --include ${DATA_TYPES} --filename "${TAXON_NAME}_ncbi.zip" || { echo "Download failed"; exit 1; }

    if [ -f "${TAXON_NAME}_ncbi.zip" ]; then
      echo "Unzipping downloaded file..."
      unzip "${TAXON_NAME}_ncbi.zip" -d "${TAXON_NAME}_ncbi" || { echo "Unzipping failed"; exit 1; }
      
      echo "Hydrating dataset..."
      datasets rehydrate --directory "${TAXON_NAME}_ncbi" || { echo "Rehydration failed"; exit 1; }
    else
      echo "Error: The file ${TAXON_NAME}_ncbi.zip was not found. The download might have failed."
      exit 1
    fi
  else
    echo "Data for ${TAXON_NAME} already downloaded."
  fi
}

# Function to check if an accession is from RefSeq
is_refseq_accession() {
  local accession="$1"
  [[ $accession =~ ^(AC_|NC_|NG_|NM_|NP_|NR_|NT_|NW_|XM_|XP_|XR_|ZP_) ]]
}

# Function to create a mapping of accessions to species names
create_species_mapping() {
  local report_file="$1"
  declare -gA species_map

  while IFS= read -r line; do
    accession=$(echo "$line" | jq -r '.accession')
    species=$(echo "$line" | jq -r '.organism.organismName' | tr ' ' '_')
    species_map["$accession"]="$species"
  done < "$report_file"
}

# Function to process individual mRNA files
process_mrna_file() {
  local dir="$1"
  local accession="$2"
  local mrna_file="$3"
  local temp_file="$4"

  species="${species_map[$accession]}"

  if [ -f "$mrna_file" ]; then
    while read -r line; do
      if [[ $line == ">"* ]]; then
        accession=$(echo "$line" | cut -d' ' -f1 | sed 's/>//')
        header=">${species}_${accession} ${line#*>}"
        if [ "$REFSEQ_ONLY" = false ] || is_refseq_accession "$accession"; then
          echo "${header}" >> "$temp_file"
        fi
        if [ "$REFSEQ_ONLY" = true ]; then
          echo "${header}" >> "$temp_file"
        fi
      else
        if [ "$REFSEQ_ONLY" = false ] || is_refseq_accession "$accession"; then
          echo "$line" >> "$temp_file"
        fi
        if [ "$REFSEQ_ONLY" = true ]; then
          echo "$line" >> "$temp_file"
        fi
      fi
    done < "$mrna_file"
  fi
}

# Function to process individual protein files
process_protein_file() {
  local dir="$1"
  local accession="$2"
  local protein_file="$3"
  local temp_file="$4"

  species="${species_map[$accession]}"

  if [ -f "$protein_file" ]; then
    while read -r line; do
      if [[ $line == ">"* ]]; then
        accession=$(echo "$line" | cut -d' ' -f1 | sed 's/>//')
        header=">${species}_${accession} ${line#*>}"
        if [ "$REFSEQ_ONLY" = false ] || is_refseq_accession "$accession"; then
          echo "${header}" >> "$temp_file"
        fi
        if [ "$REFSEQ_ONLY" = true ]; then
          echo "${header}" >> "$temp_file"
        fi
      else
        if [ "$REFSEQ_ONLY" = false ] || is_refseq_accession "$accession"; then
          echo "$line" >> "$temp_file"
        fi
        if [ "$REFSEQ_ONLY" = true ]; then
          echo "$line" >> "$temp_file"
        fi
      fi
    done < "$protein_file"
  fi
}

# Function to process files and concatenate sequences
process_files() {
  echo "Initializing output files for concatenated sequences..."
  > concatenated_mrna_refseq.fasta
  > concatenated_protein_refseq.fasta

  if [ "$REFSEQ_ONLY" = true ]; then
    > concatenated_mrna.fasta
    > concatenated_protein.fasta
  fi

  export -f process_mrna_file
  export -f process_protein_file
  export -f is_refseq_accession
  export -f process_files_parallel

  # Create the species mapping
  create_species_mapping "${TAXON_NAME}_ncbi/ncbi_dataset/data/assembly_data_report.jsonl"

  temp_dir=$(mktemp -d)
  job_count=0
  find "${TAXON_NAME}_ncbi/ncbi_dataset/data/" -type d -mindepth 1 -maxdepth 1 | while IFS= read -r dir; do
    temp_file_mrna="${temp_dir}/$(basename "$dir")_mrna.fasta"
    temp_file_protein="${temp_dir}/$(basename "$dir")_protein.fasta"
    process_files_parallel "$dir" "$temp_file_mrna" "$temp_file_protein" &
    job_count=$((job_count + 1))
    if (( job_count >= NUM_CORES )); then
      wait
      job_count=0
    fi
  done
  wait

  echo "Concatenating temporary files..."
  cat "${temp_dir}"/*_mrna.fasta >> concatenated_mrna_refseq.fasta
  cat "${temp_dir}"/*_protein.fasta >> concatenated_protein_refseq.fasta

  if [ "$REFSEQ_ONLY" = true ]; then
    cat "${temp_dir}"/*_mrna.fasta >> concatenated_mrna.fasta
    cat "${temp_dir}"/*_protein.fasta >> concatenated_protein.fasta
  fi

  rm -rf "$temp_dir"
}

# Function to process files in parallel
process_files_parallel() {
  local dir="$1"
  local temp_file_mrna="$2"
  local temp_file_protein="$3"
  local accession=$(basename "$dir")
  local mrna_file=$(find "$dir" -name "*rna.fna")
  local protein_file=$(find "$dir" -name "*protein.faa")
  process_mrna_file "$dir" "$accession" "$mrna_file" "$temp_file_mrna"
  process_protein_file "$dir" "$accession" "$protein_file" "$temp_file_protein"
}

# Function to create BLAST databases
create_blast_databases() {
  if [ -s concatenated_mrna_refseq.fasta ]; then
    echo "Creating BLAST database for mRNA..."
    makeblastdb -in concatenated_mrna_refseq.fasta -dbtype nucl -out mrna_blast_db_refseq -parse_seqids || { echo "BLAST database creation for mRNA failed"; exit 1; }
  fi

  if [ -s concatenated_protein_refseq.fasta ]; then
    echo "Creating BLAST database for proteins..."
    makeblastdb -in concatenated_protein_refseq.fasta -dbtype prot -out protein_blast_db_refseq -parse_seqids || { echo "BLAST database creation for proteins failed"; exit 1; }
  fi

  if [ "$REFSEQ_ONLY" = true ]; then
    if [ -s concatenated_mrna.fasta ]; then
      echo "Creating BLAST database for unfiltered mRNA..."
      makeblastdb -in concatenated_mrna.fasta -dbtype nucl -out mrna_blast_db -parse_seqids || { echo "BLAST database creation for unfiltered mRNA failed"; exit 1; }
    fi

    if [ -s concatenated_protein.fasta ]; then
      echo "Creating BLAST database for unfiltered proteins..."
      makeblastdb -in concatenated_protein.fasta -dbtype prot -out protein_blast_db -parse_seqids || { echo "BLAST database creation for unfiltered proteins failed"; exit 1; }
    fi
  fi
}

# Main script execution
main() {
  download_data
  process_files
  create_blast_databases

  echo "Download and processing complete. The NCBI data for ${TAXON_NAME} with data types: ${DATA_TYPES}, including concatenated mRNA and protein files and BLAST databases, is available."
}

# Run the main function
main
