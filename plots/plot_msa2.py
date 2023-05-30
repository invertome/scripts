import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from matplotlib.colors import ListedColormap
import json

def parse_positions(positions_str):
    positions = []
    for pos_str in positions_str.split(','):
        if '-' in pos_str:  # If there's a dash, it represents a range
            start, end = map(int, pos_str.split('-'))  # Get the start and end of the range
            positions.extend(range(start, end + 1))  # Add all positions in the range to the list
        elif pos_str:  # If there's no dash and the string is not empty, it's a single position
            positions.append(int(pos_str))  # Add the position to the list
    return positions  # Return the list of positions

def parse_range(range_str):
    try:
        start, end = map(int, range_str.split('-'))
    except ValueError:
        print("Invalid range string. Range must be two integers separated by a dash.")
        return None, None

    if start > end:
        print("Invalid range specified. Please make sure your range consists of two numbers, with the first being smaller or equal to the second.")
        return None, None

    return start, end

def create_color_map_dict(color_scheme_path):
    # Load the color scheme from the JSON file
    with open(color_scheme_path, 'r') as f:
        color_scheme = json.load(f)
    color_map_dict = color_scheme["colors"]
    return color_map_dict

def plot_msa(input_file, output_folder, plot_range, highlight_positions, color_map_dict):
    # Read the alignment from the input file
    alignment = AlignIO.read(input_file, "fasta")
    # Select only the sequences in the specified range
    alignment = alignment[:, plot_range[0]-1:plot_range[1]]

    # Convert the alignment to a matrix for easier handling
    matrix = np.array([list(rec) for rec in alignment], dtype='U1')

    # Create a figure and axis for the plot
    fig, ax = plt.subplots(figsize=(10, len(alignment)*0.3))

    # Plot the grid
    for pos in range(matrix.shape[1]):
        for seq in range(matrix.shape[0]):
            ax.text(pos+0.5, seq+0.5, matrix[seq, pos], color=color_map_dict.get(matrix[seq, pos], 'black'), ha='center', va='center')

    # Highlight the specified positions
    for pos in highlight_positions:
        if plot_range[0] <= pos <= plot_range[1]: 
            ax.axvspan(xmin=pos - plot_range[0], xmax=pos - plot_range[0] + 1, color='red', alpha=0.3)

    # Set the labels for the axes
    ax.set_xticks(np.arange(matrix.shape[1])+0.5)
    ax.set_xticklabels(np.arange(plot_range[0], plot_range[1] + 1))
    ax.set_yticks(np.arange(matrix.shape[0])+0.5)
    ax.set_yticklabels([rec.id for rec in alignment], rotation=0)
    ax.set_title(os.path.basename(input_file))

    # Remove axes
    ax.axis('off')
    
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
                        help='Positions for highlighting, comma separated. Ranges can be specified with a dash.')
    parser.add_argument('-r', '--range', type=str, required=True,
                        help='Range of positions to plot, two integers separated by a dash.')
    parser.add_argument('-c', '--color_scheme', type=str, required=True,
                        help='Path to the color scheme JSON file.')
    args = parser.parse_args()

    # Parse the positions for highlighting and the range for plotting
    highlight_positions = parse_positions(args.positions) if args.positions else []
    start, end = parse_range(args.range)

    # Check if the range is valid
    if start is None or end is None:
        return
    plot_range = (start, end)

    # Get the color map dictionary from the color scheme JSON file
    color_map_dict = create_color_map_dict(args.color_scheme)

    # Get the base name of the file (without extension)
    base_name = os.path.splitext(os.path.basename(args.input_file))[0]
    # Create the output folder based on the base name of the file
    output_folder = os.path.join(os.getcwd(), base_name)
    os.makedirs(output_folder, exist_ok=True)

    plot_msa(args.input_file, output_folder, plot_range, highlight_positions, color_map_dict)

if __name__ == '__main__':
    main()
