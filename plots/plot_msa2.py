import argparse
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from matplotlib.patches import Rectangle
from matplotlib import rcParams

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

def plot_msa(input_file, output_folder, plot_range, highlight_positions, color_map_dict):
    # Read the alignment from the input file
    alignment = AlignIO.read(input_file, "fasta")
    # Select only the sequences in the specified range
    alignment = alignment[:, plot_range[0]-1:plot_range[1]]

    # Convert the alignment to a matrix for easier handling
    matrix = np.array([list(rec) for rec in alignment], dtype='U1')

    # Create a figure and axis for the plot
    fig, ax = plt.subplots(figsize=(len(alignment[0]), len(alignment)*0.3))

    # Plot the grid
    for pos in range(matrix.shape[1]):
        for seq in range(matrix.shape[0]):
            ax.text(pos+0.5, seq+0.5, matrix[seq, pos], color='black', ha='center', va='center', fontsize=8)
            ax.add_patch(Rectangle((pos, seq), 1, 1, color=color_map_dict.get(matrix[seq, pos], 'white')))

    # Highlight the specified positions
    for pos in highlight_positions:
        if plot_range[0] <= pos <= plot_range[1]: 
            ax.add_patch(Rectangle((pos - plot_range[0], -0.5), 1, matrix.shape[0], color='red', alpha=0.3))

    # Set the labels for the axes
    ax.set_xticks(np.arange(matrix.shape[1])+0.5)
    ax.set_xticklabels(np.arange(plot_range[0], plot_range[1] + 1), fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0])+0.5)
    ax.set_yticklabels([rec.id for rec in alignment], rotation=0, fontsize=8)
    ax.set_title(os.path.basename(input_file), fontweight='bold', pad=20)
    
    # Adjust layout
    fig.tight_layout()

    # Save the plot as PDF and PNG
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + ".pdf"))
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + ".png"))
    plt.close()  # Close the plot

def main():
    # Define the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True,
                        help='Input fasta file.')
    parser.add_argument('-p', '--positions', type=str, default='',
                        help='Positions to highlight. Can be a single number, a range (e.g., "1-3"), or multiple numbers/ranges separated by commas (e.g., "1,3,5-7").')
    parser.add_argument('-r', '--range', type=str, required=True,
                        help='Range of the alignment to plot, e.g., "1-100".')
    parser.add_argument('-c', '--color_map', type=str, required=True,
                        help='JSON file with color map for amino acids.')
    parser.add_argument('-o', '--output_folder', type=str, default='.',
                        help='Folder to save the output plot.')

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
