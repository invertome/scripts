import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Function to search for HMM motifs using the hmmsearch command
def search_hmm_profile(hmm_file, input_fasta, output_table, e_value=None, threads=1):
    command = f"hmmsearch --tblout {output_table} "
    
    if e_value is not None:
        command += f"--domE {e_value} "
    
    command += f"--cpu {threads} {hmm_file} {input_fasta}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

    with open(output_table, "r") as f:
        lines = [line for line in f if not line.startswith("#")]
        return len(lines)

# Function to filter the results by E-value
def filter_results_by_evalue(input_table, output_table, e_value):
    with open(input_table, "r") as fin, open(output_table, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                fields = line.strip().split()
                hit_e_value = float(fields[4])
                if hit_e_value <= e_value:
                    fout.write(line)

# Function to parse the species tab file
def parse_species_tab(species_tab_file):
    species_df = pd.read_csv(species_tab_file, sep="\t", header=None)
    species_df.columns = ["id", "species_id", "species_name", "other1", "other2", "other3", "other4"]
    species_dict = species_df.set_index("species_id")["species_name"].to_dict()
    return species_dict

# Function to read the results table and add species information
def read_results_table(output_table, species_dict=None):
    with open(output_table, "r") as f:
        # Count the number of initial comment lines to skip them when reading the file
        skip_rows = 0
        line = f.readline()
        while line.startswith("#"):
            skip_rows += 1
            line = f.readline()
            
    print(f"Output table: {output_table}")
    print(f"Number of lines to skip: {skip_rows}")

    # Read the file manually to handle varying number of columns
    data = []
    with open(output_table, "r") as f:
        for _ in range(skip_rows):
            f.readline()
        for line in f:
            columns = line.strip().split(maxsplit=18)  # Split into the first 18 columns
            if len(columns) > 18:
                columns[18] = " ".join(columns[18:])  # Join remaining columns into the description
                columns = columns[:19]
            data.append(columns)

    column_names = ["target_name", "target_accession", "query_name", "query_accession", "E-value", "score", "bias", "domain_E-value", "domain_score", "domain_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description"]

    results_df = pd.DataFrame(data, columns=column_names)

    print("Loaded dataframe:")
    print(results_df.head())
    print(results_df.columns)

    # Ensure E-value column is numeric
    results_df["E-value"] = pd.to_numeric(results_df["E-value"], errors='coerce')

    if species_dict and not results_df.empty:
        results_df["species_id"] = results_df["target_name"].apply(lambda x: x.split(":")[-1].split("_")[0])
        results_df["species_name"] = results_df["species_id"].apply(lambda x: species_dict.get(x, "Unknown"))
        results_df["description_with_species"] = results_df["species_name"] + "_" + results_df["description"]

    return results_df

# Function to plot the E-value distribution of hits
def plot_evalue_distribution(results_df):
    plt.figure()
    plt.hist(np.log10(results_df["E-value"].dropna()), bins=30)
    plt.xlabel("log10(E-value)")
    plt.ylabel("Frequency")
    plt.title("E-value distribution of hits")
    plt.savefig("evalue_distribution.png")
    plt.close()

# Function to plot the taxonomic distribution of hits
def plot_taxonomic_distribution(results_df):
    if "species_name" not in results_df.columns:
        print("species_name column not found in the dataframe.")
        return
    plt.figure()
    species_counts = results_df["species_name"].value_counts()
    species_counts.plot(kind="barh")
    plt.xlabel("Number of hits")
    plt.ylabel("Species")
    plt.title("Taxonomic distribution of hits")
    plt.savefig("taxonomic_distribution.png")
    plt.close()

# Function to plot the E-value threshold vs. number of hits
def plot_threshold_vs_hits(thresholds, num_hits):
    plt.figure()
    plt.plot(thresholds, num_hits, marker="o")
    plt.xlabel("E-value threshold")
    plt.ylabel("Number of hits")
    plt.title("E-value threshold vs. number of hits")
    plt.savefig("threshold_vs_hits.png")
    plt.close()

def main(args):
    if args.species_tab_file:
        species_dict = parse_species_tab(args.species_tab_file)
    else:
        species_dict = None

    if args.plot:
        output_table_initial = f"{args.output}_initial.txt"

        # Perform an initial permissive search to determine the range of E-values
        search_hmm_profile(args.hmm_file, args.input, output_table_initial, e_value=10, threads=args.threads)
        initial_hits_df = read_results_table(output_table_initial, species_dict)
        
        # Determine the min and max E-values from the initial search
        min_evalue = initial_hits_df["E-value"].min() if not initial_hits_df.empty else 1e-50
        max_evalue = initial_hits_df["E-value"].max() if not initial_hits_df.empty else 1e-1
        
        # Generate E-values for the searches
        e_values = np.logspace(np.log10(min_evalue), np.log10(max_evalue), num=args.searches)[::-1]
        print(f"Generated E-values: {e_values}")
        num_hits = []

        output_table_most_permissive = f"{args.output}_most_permissive.txt"

        if not args.use_existing_results:
            # Perform search with the most permissive E-value
            search_hmm_profile(args.hmm_file, args.input, output_table_most_permissive, e_value=e_values[0], threads=args.threads)

        for e_value in e_values:
            output_table_filtered = f"{args.output}_evalue_{e_value:.1e}.txt"
            filter_results_by_evalue(output_table_most_permissive, output_table_filtered, e_value)
            hits_df = read_results_table(output_table_filtered, species_dict)
            num_hits.append(hits_df.shape[0])
            print(f"E-value: {e_value}, Number of hits: {hits_df.shape[0]}")

        plot_threshold_vs_hits(e_values, num_hits)
        plot_evalue_distribution(hits_df)
        plot_taxonomic_distribution(hits_df)

    else:
        if not args.use_existing_results:
            search_hmm_profile(args.hmm_file, args.input, args.output, threads=args.threads)

        if species_dict:
            hits_df = read_results_table(args.output, species_dict)
            plot_evalue_distribution(hits_df)
            plot_taxonomic_distribution(hits_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Search for HMM motifs in a FASTA file.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output results file")
    parser.add_argument("-m", "--hmm_file", type=str, required=True, help="HMM profile file")
    parser.add_argument("-p", "--plot", action="store_true", help="Plot threshold vs. number of hits")
    parser.add_argument("-n", "--searches", type=int, default=10, help="Number of E-values to search")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads/CPUs to use")
    parser.add_argument("-x", "--species_tab_file", type=str, help="Species tab file for taxonomic information")
    parser.add_argument("-k", "--use_existing_results", action="store_true", help="Use existing results file")
    args = parser.parse_args()
    main(args)
