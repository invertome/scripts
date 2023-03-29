# USAGE: python correct_newick.py -i input_file.tre -o output_file.tre

import re
import argparse
from ete3 import Tree

def correct_newick_tree(input_file, output_file):
    # Read the input Newick tree file
    with open(input_file, 'r') as f:
        tree_string = f.read()

    # Find and correct redundant double-brackets
    corrected_tree_string = re.sub(r'\(\(([^\(\)]+)\)\)', r'(\1)', tree_string)

    # Check if any more redundant double-brackets are found and correct them
    while re.search(r'\(\(([^\(\)]+)\)\)', corrected_tree_string):
        corrected_tree_string = re.sub(r'\(\(([^\(\)]+)\)\)', r'(\1)', corrected_tree_string)

    try:
        # Parse the original and corrected Newick trees
        original_tree = Tree(tree_string)
        corrected_tree = Tree(corrected_tree_string)

        # Check if the topology of the corrected tree is the same as the original tree
        if original_tree.robinson_foulds(corrected_tree, unrooted_trees=True)[0] == 0:
            # Write the corrected tree to the output file
            with open(output_file, 'w') as f:
                f.write(corrected_tree_string)
            print(f"Corrected Newick tree written to {output_file}")
        else:
            # Warn the user if the topology differs
            print("Warning: Corrected tree topology differs from the original tree. Manual inspection is required.")

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
