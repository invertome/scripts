import argparse
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
from Bio import AlignIO
import logomaker as lm
from collections import Counter
import math

def parse_positions(positions_str):
    """Parse a string of positions into a list of integers."""
    positions = []
    for pos_str in positions_str.split(','):
        if '-' in pos_str:
            start, end = map(int, pos_str.split('-'))
            positions.extend(range(start, end + 1))
        elif pos_str:
            positions.append(int(pos_str))
    return positions

def load_color_map_from_json(json_file):
    """Load color map from a JSON file with error handling."""
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Color map file {json_file} not found.")
    except json.JSONDecodeError:
        raise ValueError(f"Color map file {json_file} is not a valid JSON.")
    color_map_dict = data['colors']
    color_map_dict['-'] = 'white'
    return color_map_dict

def compute_counts(matrix):
    """Compute amino acid counts for sequence logo."""
    alphabet = 'ACDEFGHIKLMNPQRSTVWY-'
    counts_df = pd.DataFrame(0, index=range(matrix.shape[1]), columns=list(alphabet))
    for pos in range(matrix.shape[1]):
        counter = Counter(matrix[:, pos])
        for aa, count in counter.items():
            if aa in counts_df.columns:
                counts_df.loc[pos, aa] = count
    return counts_df

def plot_msa(input_file, output_folder, plot_range, highlight_positions, color_map_dict, wrap=100):
    """Plot multiple sequence alignment with sequence logo."""
    # Load alignment and compute original ranges
    alignment = AlignIO.read(input_file, "fasta")
    original_ranges = []
    for seq in alignment:
        seq_str = str(seq.seq)
        non_gap_indices = [i + 1 for i, char in enumerate(seq_str) if char != '-']
        if non_gap_indices:
            original_range = (non_gap_indices[0] - non_gap_indices.count(1) + 1,
                              non_gap_indices[-1] - non_gap_indices.count(1) + 1)
        else:
            original_range = (None, None)
        seq_id = seq.id.split('/')[0] if '/' in seq.id else seq.id
        original_ranges.append((seq_id.replace('_', ' '), original_range))
    
    # Slice alignment to the specified range
    alignment = alignment[:, plot_range[0] - 1:plot_range[1]]
    block_size = wrap
    total_positions = len(alignment[0])
    number_of_blocks = math.ceil(total_positions / block_size)
    n_seqs = len(alignment)
    
    # Allocate more vertical space for the logo plot
    logo_height = 3  # Increase logo plot height (in inches) to avoid squishing
    fig_width = block_size * 0.3  # Each column is 0.3 inches wide
    fig_height = (n_seqs * 0.3 + logo_height) * number_of_blocks
    fig, axs = plt.subplots(2 * number_of_blocks, 1,
                            figsize=(fig_width, fig_height),
                            gridspec_kw={'height_ratios': [n_seqs, logo_height] * number_of_blocks, 'hspace': 0.5})
    fig.subplots_adjust(top=0.95)  # Reduce blank space above first block
    if number_of_blocks == 1:
        axs = np.array([axs]).flatten()
    
    # Plot each block
    for k in range(number_of_blocks):
        start_idx = k * block_size
        end_idx = min(start_idx + block_size, total_positions)
        matrix = np.array([list(str(rec.seq)[start_idx:end_idx]) for rec in alignment], dtype='U1')
        counts_df = compute_counts(matrix)
        
        # Alignment subplot and corresponding logo subplot
        ax = axs[2 * k]
        ax_logo = axs[2 * k + 1]
        
        # Convert matrix to numerical for imshow
        alphabet = list('ACDEFGHIKLMNPQRSTVWY-')
        aa_to_index = {aa: i for i, aa in enumerate(alphabet)}
        numerical_matrix = np.array([[aa_to_index.get(aa, -1) for aa in row] for row in matrix])
        colors = [color_map_dict.get(aa, 'white') for aa in alphabet]
        cmap = matplotlib.colors.ListedColormap(colors)
        # Set extent so that each cell is 1x1 and centers are at (pos+0.5, seq+0.5)
        ax.imshow(numerical_matrix, aspect='auto', cmap=cmap, interpolation='none',
                  origin='upper', extent=[0, matrix.shape[1], matrix.shape[0], 0])
        
        # Add text overlay with corrected positioning
        for pos in range(matrix.shape[1]):
            actual_pos = plot_range[0] + start_idx + pos
            for seq in range(matrix.shape[0]):
                if actual_pos in highlight_positions:
                    ax.text(pos + 0.5, seq + 0.5, matrix[seq, pos],
                            color='black', ha='center', va='center', fontsize=10, fontweight='bold')
                else:
                    ax.text(pos + 0.5, seq + 0.5, matrix[seq, pos],
                            color='black', ha='center', va='center', fontsize=8)
        
        # Set axis limits and ticks; note the reversed y-axis to match the extent
        ax.set_xlim(0, matrix.shape[1])
        ax.set_ylim(matrix.shape[0], 0)
        ax.set_xticks(np.arange(matrix.shape[1]) + 0.5)
        ax.set_xticklabels(range(plot_range[0] + start_idx, plot_range[0] + start_idx + matrix.shape[1]),
                           rotation=90, fontsize=6)
        
        # Set y-axis ticks and labels (sequence names)
        ax.set_yticks(np.arange(matrix.shape[0]) + 0.5)
        ax.set_yticklabels([f"{seq_id} " for seq_id, rng in original_ranges],
                           fontsize=8, va='center')
        ax.tick_params(axis='both', which='both', length=0)
        
        # Plot sequence logo
        counts_df = counts_df.divide(counts_df.sum(axis=1), axis=0)
        counts_df.index = range(plot_range[0] + start_idx, plot_range[0] + start_idx + matrix.shape[1])
        lm.Logo(counts_df, ax=ax_logo, color_scheme=color_map_dict)
    
    # Add title and save plot
    fig.suptitle(f'Multiple Sequence Alignment of {os.path.basename(input_file)}',
                 fontsize=12, fontweight='bold')
    output_folder = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0])
    os.makedirs(output_folder, exist_ok=True)
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + ".pdf"), bbox_inches='tight')
    plt.savefig(os.path.join(output_folder, os.path.basename(input_file) + ".png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save color scheme
    with open(os.path.join(output_folder, "color_scheme.txt"), "w") as file:
        for key, value in color_map_dict.items():
            file.write(f"{key} : {value}\n")

def main():
    """Parse arguments and run the plotting function."""
    parser = argparse.ArgumentParser(description='This script generates a plot for a multiple sequence alignment. It provides options to highlight specific positions and to use a custom color map.')
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='Path to the input file containing the multiple sequence alignment in FASTA format.')
    parser.add_argument('-p', '--positions', type=str, required=False, default='',
                        help='Positions in the alignment to highlight. This can be a single number, a range (e.g., "1-3"), or multiple numbers/ranges separated by commas (e.g., "1,3,5-7").')
    parser.add_argument('-r', '--range', type=str, required=True,
                        help='Range of positions in the alignment to plot. This should be in the format "start-end", e.g., "1-100".')
    parser.add_argument('-c', '--color_map', type=str, required=True,
                        help='Path to a JSON file containing a color map for the amino acids.')
    parser.add_argument('-o', '--output_folder', type=str, default='.',
                        help='Path to the folder where the output plot should be saved.')
    parser.add_argument('-w', '--wrap', type=int, default=100,
                        help='Number of positions per block')
    args = parser.parse_args()
    
    # Error handling for input file and range
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file {args.input_file} not found.")
    range_start, range_end = map(int, args.range.split('-'))
    if range_start > range_end:
        raise ValueError("Range start must be less than or equal to range end.")
    
    highlight_positions = parse_positions(args.positions)
    color_map_dict = load_color_map_from_json(args.color_map)
    plot_range = (range_start, range_end)
    plot_msa(args.input_file, args.output_folder, plot_range, highlight_positions, color_map_dict, args.wrap)

if __name__ == "__main__":
    main()
