import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt

def search_hmm_profile(hmm_file, input_fasta, output_table, e_value=None, threads=1):
    command = f"hmmsearch --tblout {output_table} "
    
    if e_value is not None:
        command += f"--incE {e_value} "
    
    command += f"--cpu {threads} {hmm_file} {input_fasta}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

    with open(output_table, "r") as f:
        lines = [line for line in f if not line.startswith("#")]
        return len(lines)

def plot_threshold_vs_hits(thresholds, num_hits):
    plt.plot(thresholds, num_hits, marker="o")
    plt.xlabel("E-value threshold")
    plt.ylabel("Number of hits")
    plt.xscale("log")
    plt.title("E-value threshold vs. number of hits")
    plt.savefig("threshold_vs_hits.png")

def main(args):
    if not args.plot_thresholds:
        search_hmm_profile(args.hmm_file, args.input_fasta, args.output_table, threads=args.threads)
    else:
        e_value_thresholds = np.logspace(-args.start, -args.end, args.searches)
        num_hits = []
        for e_value in e_value_thresholds:
            output_table = f"results_e{e_value:.3e}.out"
            hits = search_hmm_profile(
                args.hmm_file, args.input_fasta, output_table, e_value=e_value, threads=args.threads
            )
            num_hits.append(hits)
        plot_threshold_vs_hits(e_value_thresholds, num_hits)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Search HMM profile in a FASTA file.")
    parser.add_argument("-i", "--input_fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output_table", required=True, help="Output table file")
    parser.add_argument("-m", "--hmm_file", required=True, help="HMM profile file")
    parser.add_argument(
        "-p",
        "--plot_thresholds",
        action="store_true",
        help="Plot E-value threshold vs. number of hits",
    )
    parser.add_argument(
        "-s", "--start", type=int, default=1, help="Start positive exponent for E-value range"
    )
    parser.add_argument(
        "-e", "--end", type=int, default=10, help="End positive exponent for E-value range"
    )
    parser.add_argument(
        "-n", "--searches", type=int, default=10, help="Number of E-values to search"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads/CPUs to use"
    )
    args = parser.parse_args()
    main(args)
