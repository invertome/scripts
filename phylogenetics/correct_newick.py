# USAGE: python correct_newick.py -i input_file.tre -o output_file.tre

import argparse
from io import StringIO
from Bio import Phylo

def collapse_single_child_nodes(tree):
    for clade in tree.find_clades():
        if len(clade.clades) == 1:
            clade.clades[0].branch_length += clade.branch_length
            clade.clades = clade.clades[0].clades
            clade.name = clade.clades[0].name
            clade.confidence = clade.clades[0].confidence
    return tree

def correct_newick_tree(input_file, output_file):
    try:
        # Parse the original Newick tree
        original_tree = Phylo.read(input_file, 'newick')

        # Collapse nodes with a single child
        corrected_tree = collapse_single_child_nodes(original_tree)

        # Write the corrected tree to the output file
        Phylo.write(corrected_tree, output_file, 'newick')
        print(f"Corrected Newick tree written to {output_file}")

    except Exception as e:
        # Print error message if there's an issue while parsing the trees
        print(f"Error while parsing trees: {str(e)}")
        print("Manual inspection of the input file is required.")

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description='Correct Newick tree syntax issues')
parser.add_argument('-i', '--input', type=str, required=True, help='Path to input Newick tree file')
parser.add_argument('-o', '--output', type=str, required=True, help='Path to corrected output Newick tree file')
args = parser.parse_args()

# Call the correct_newick_tree function with the specified input and output file paths
correct_newick_tree(args.input, args.output)

