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

    try:
        results_df = pd.read_csv(output_table, sep="\s+", skiprows=skip_rows, header=None)
        column_names = ["target_name", "target_accession", "query_name", "query_accession", "E-value", "score", "bias", "domain_E-value", "domain_score", "domain_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description"]
        results_df.columns = column_names
    except pd.errors.EmptyDataError:
        column_names = ["target_name", "target_accession", "query_name", "query_accession", "E-value", "score", "bias", "domain_E-value", "domain_score", "domain_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description"]
        results_df = pd.DataFrame(columns=column_names)

    print("Loaded dataframe:")
    print(results_df.head())
    print(results_df.columns)


    if species_dict and not results_df.empty:
        results_df["species_id"] = results_df["target_name"].apply(lambda x: x.split(":")[-1])
        results_df["species_name"] = results_df["species_id"].apply(lambda x: species_dict.get(x, "Unknown"))
        results_df["description_with_species"] = results_df["species_name"] + "_" + results_df["description"]

    return results_df


# Function to plot the E-value distribution of hits
def plot_evalue_distribution(results_df):
    plt.figure()
    plt.hist(np.log10(results_df["E-value"]), bins=30)
    plt.xlabel("log10(E-value)")
    plt.ylabel("Frequency")
    plt.title("E-value distribution of hits")
    plt.savefig("evalue_distribution.png")

# Function to plot the taxonomic distribution of hits
def plot_taxonomic_distribution(results_df):
    plt.figure()
    species_counts = results_df["species_name"].value_counts()
    species_counts.plot(kind="barh")
    plt.xlabel("Number of hits")
    plt.ylabel("Species")
    plt.title("Taxonomic distribution of hits")
    plt.savefig("taxonomic_distribution.png")

# Function to plot the E-value threshold vs. number of hits
def plot_threshold_vs_hits(thresholds, num_hits):
    plt.plot(thresholds, num_hits, marker="o")
    plt.xlabel("E-value threshold")
    plt.ylabel("Number of hits")
    plt.title("E-value threshold vs. number of hits")
    plt.savefig("threshold_vs_hits.png")

def main(args):
    if args.species_tab_file:
        species_dict = parse_species_tab(args.species_tab_file)
    else:
        species_dict = None

    if args.plot:
        e_values = np.logspace(-args.end, -args.start, num=args.searches)[::-1]
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
    parser.add_argument("-s", "--start", type=float, default=1, help="Start log10(E-value) for plot")
    parser.add_argument("-e", "--end", type=float, default=50, help="End log10(E-value) for plot")
    parser.add_argument("-n", "--searches", type=int, default=10, help="Number of E-values to search")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads/CPUs to use")
    parser.add_argument("-x", "--species_tab_file", type=str, help="Species tab file for taxonomic information")
    parser.add_argument("-k", "--use_existing_results", action="store_true", help="Use existing results file")
    args = parser.parse_args()
    main(args)
