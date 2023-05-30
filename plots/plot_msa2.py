import argparse
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import AlignIO
from matplotlib.patches import Rectangle
from matplotlib import rcParams
import logomaker as lm
from collections import Counter

def parse_positions(positions_str):
    positions = []
    for pos_str in positions_str.split(','):
        if '-' in pos_str:  # If there's a dash, it represents a range
            start, end = map(int, pos_str.split('-'))  # Get the start and end of the range
            positions.extend(range(start, end + 1))  # Add all positions in the range to the list
        elif pos_str:  # If there's no dash and the string is not empty, it's a single position
            positions.append(int(pos_str))  # Add the position to the list
    return positions  # Return the list of positions

def load_color_map_from_json(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)
    color_map_dict = data['colors']
    return color_map_dict

def compute_counts(matrix):
    alphabet = 'ACDEFGHIKLMNPQRSTVWY-'
    counts_df = pd.DataFrame(0, index=range(matrix.shape[1]), columns=list(alphabet))
    for pos in range(matrix.shape[1]):
        counter = Counter(matrix[:, pos])
        for aa, count in counter.items():
            if aa in counts_df.columns:
                counts_df.loc[pos, aa] = count
    return counts_df

def plot_msa(input_file, output_folder, plot_range, highlight_positions, color_map_dict):
    # Read the alignment from the input file
    alignment = AlignIO.read(input_file, "fasta")
    # Select only the sequences in the specified range
    alignment = alignment[:, plot_range[0]-1:plot_range[1]]

    # Convert the alignment to a matrix for easier handling
    matrix = np.array([list(rec) for rec in alignment], dtype='U1')

    # Compute counts for logo plot
    counts_df = compute_counts(matrix)

    # Create a figure and axis for the plot
    fig, axs = plt.subplots(2, 1, figsize=(len(alignment[0])*0.3, len(alignment)*0.3), gridspec_kw={'height_ratios': [3, 1]})
    ax, ax_logo = axs

    # Plot the grid
    for pos in range(matrix.shape[1]):
        for seq in range(matrix.shape[0]):
            if pos + plot_range[0] in highlight_positions:  # if position needs to be highlighted
                ax.text(pos+0.5, matrix.shape[0]-seq-0.5, matrix[seq, pos], color='black', ha='center', va='center', fontsize=10, fontweight='bold')
            else:
                ax.text(pos+0.5, matrix.shape[0]-seq-0.5, matrix[seq, pos], color='black', ha='center', va='center', fontsize=8)
            ax.add_patch(Rectangle((pos, matrix.shape[0]-seq-1), 1, 1, edgecolor='gray', facecolor=color_map_dict.get(matrix[seq, pos], 'white'), lw=0))

    # Set the axis limits
    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(0, matrix.shape[0])
    
    # Add x and y axis labels
    ax.set_xticks(np.arange(matrix.shape[1])+0.5)
    ax.set_xticklabels(range(plot_range[0], plot_range[1]+1), rotation=90, fontsize=6)
    ax.set_yticks(np.arange(matrix.shape[0])+0.5)
    ax.set_yticklabels([rec.id for rec in reversed(alignment)], fontsize=8, va='center')

    # Turn off the tick marks
    ax.tick_params(axis=u'both', which=u'both',length=0)

    # Set the title
    ax.set_title(f'Multiple Sequence Alignment of {os.path.basename(input_file)}', fontsize=12, fontweight='bold', pad=20)

    # Convert counts to percentages for sequence logo
    counts_df = counts_df.divide(counts_df.sum(axis=1), axis=0)
    
    # Make a sequence logo
    counts_df.index = range(plot_range[0], plot_range[1]+1)
    lm.Logo(counts_df, ax=ax_logo, color_scheme=color_map_dict)

    # Output to file
    output_folder = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0])
    os.makedirs(output_folder, exist_ok=True)
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + ".pdf"), bbox_inches='tight')
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + ".png"), dpi=300, bbox_inches='tight')
    plt.close()

    # Save color map used
    with open(os.path.join(output_folder, "color_scheme.txt"), "w") as file:
        for key, value in color_map_dict.items():
            file.write(f"{key} : {value}\n")


def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Plot a multiple sequence alignment.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input file with the multiple sequence alignment in FASTA format.')
    parser.add_argument('-p', '--positions', type=str, required=False, default='', help='Positions to highlight. Can be a single number, a range (e.g., "1-3"), or multiple numbers/ranges separated by commas (e.g., "1,3,5-7").')
    parser.add_argument('-r', '--range', type=str, required=True, help='Range of the alignment to plot, e.g., "1-100".')
    parser.add_argument('-c', '--color_map', type=str, required=True, help='JSON file with color map for amino acids.')
    parser.add_argument('-o', '--output_folder', type=str, default='.', help='Folder to save the output plot.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Parse the positions to highlight
    highlight_positions = parse_positions(args.positions)

    # Load the color map from the JSON file
    color_map_dict = load_color_map_from_json(args.color_map)

    # Parse the range
    range_start, range_end = map(int, args.range.split('-'))
    plot_range = range_start, range_end

    # Plot the MSA
    plot_msa(args.input_file, args.output_folder, plot_range, highlight_positions, color_map_dict)

if __name__ == "__main__":
    main()
