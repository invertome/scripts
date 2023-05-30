import argparse
import os
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns

# Function to parse the string of positions passed as argument
def parse_positions(positions_str):
    positions = []

    # Split the string by commas and iterate through the elements
    for pos_str in positions_str.split(','):
        if '-' in pos_str:  # If there's a dash, it represents a range
            start, end = map(int, pos_str.split('-'))  # Get the start and end of the range
            positions.extend(range(start, end + 1))  # Add all positions in the range to the list
        elif pos_str:  # If there's no dash and the string is not empty, it's a single position
            positions.append(int(pos_str))  # Add the position to the list

    return positions  # Return the list of positions




# Function to plot a multiple sequence alignment with a specified color map
def plot_msa(input_file, output_folder, plot_range, highlight_positions, color_map):
    # Read the alignment from the input file
    alignment = AlignIO.read(input_file, "fasta")
    # Select only the sequences in the specified range
    alignment = alignment[:, plot_range[0]-1:plot_range[1]]

    # Convert the alignment to a matrix for easier handling
    matrix = np.array([list(rec) for rec in alignment], 'S1') # changed np.character to 'S1'
    
    # Create a figure and axis for the plot
    fig, ax = plt.subplots(figsize=(10, len(alignment)*0.3))

    # Get the color map
    cmap = plt.cm.get_cmap(color_map)
    # Create a heatmap with the color map
    sns.heatmap(data=np.arange(matrix.shape[1]).reshape(1, -1), cmap=cmap, ax=ax, cbar=False)

    # Highlight the specified positions
    for pos in highlight_positions:
        if plot_range[0] <= pos <= plot_range[1]: 
            ax.axvline(x=pos - plot_range[0] + 0.5, color='red', linewidth=2)

    # Plot the grid
    for pos in range(matrix.shape[1]):
        for seq in range(matrix.shape[0]):
            ax.text(pos+1, seq, matrix[seq, pos].decode('UTF-8'), color='black', ha='center', va='center')

    # Set the labels for the axes
    ax.set_xticks(np.arange(matrix.shape[1])+0.5)
    ax.set_xticklabels(np.arange(plot_range[0], plot_range[1] + 1))
    ax.set_yticks(np.arange(matrix.shape[0])+0.5)
    ax.set_yticklabels([rec.id for rec in alignment], rotation=0)
    ax.set_title(os.path.basename(input_file))
    
    # Adjust layout
    fig.tight_layout()
    
    # Save the plot as PDF and PNG
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + "_" + color_map + ".pdf"))
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + "_" + color_map + ".png"))
    plt.close()  # Close the plot



def parse_range(range_str):
    # Split the string by the dash and map each part to an integer
    start, end = map(int, range_str.split('-'))
    return start, end

def main():
    # Define the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_files', nargs='+', required=True,
                        help='Input fasta files, space separated if multiple.')
    parser.add_argument('-p', '--positions', type=str, default='',
                        help='Positions for highlighting, comma separated. Ranges can be specified with a dash.')
    parser.add_argument('-r', '--range', type=str, required=True,
                        help='Range of positions to plot, two integers separated by a dash.')
    args = parser.parse_args()

    # Parse the positions for highlighting and the range for plotting
    highlight_positions = parse_positions(args.positions) if args.positions else []
    plot_range = parse_range(args.range)  # Use the parse_range function here

    # Check if the range is valid
    if plot_range[0] > plot_range[1]:
        print("Invalid range specified. Please make sure your range consists of two numbers, with the first being smaller or equal to the second.")
        return

   
    # Iterate through the input files
    for input_file in args.input_files:
        # Get the base name of the file (without extension)
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        # Create the output folder based on the base name of the file
        output_folder = os.path.join(os.getcwd(), base_name)
        os.makedirs(output_folder, exist_ok=True)

        # Generate plots with different color maps
        for cmap in ['viridis', 'plasma', 'inferno', 'magma']:
            plot_msa(input_file, output_folder, plot_range, highlight_positions, cmap)


if __name__ == '__main__':
    main()
