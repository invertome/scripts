# USAGE: python correct_newick.py -i input_file.tre -o output_file.tre

import re
import argparse
from ete3 import Tree

def custom_tree_parser(newick_string):
    newick_string = re.sub(r'(?<=\))(?=:\d+\.\d+)', 'missing_name', newick_string)
    return Tree(newick_string)

def correct_newick_tree(input_file, output_file):
    with open(input_file, 'r') as f:
        tree_string = f.read()

    corrected_tree_string = re.sub(r'\(\(([^\(\)]+)\)\)', r'(\1)', tree_string)

    while re.search(r'\(\(([^\(\)]+)\)\)', corrected_tree_string):
        corrected_tree_string = re.sub(r'\(\(([^\(\)]+)\)\)', r'(\1)', corrected_tree_string)

    try:
        original_tree = custom_tree_parser(tree_string)
        corrected_tree = custom_tree_parser(corrected_tree_string)

        if original_tree.robinson_foulds(corrected_tree, unrooted_trees=True)[0] == 0:
            with open(output_file, 'w') as f:
                f.write(corrected_tree_string)
            print(f"Corrected Newick tree written to {output_file}")
        else:
            print("Warning: Corrected tree topology differs from the original tree. Manual inspection is required.")

    except Exception as e:
        print(f"Error while parsing trees: {str(e)}")
        print("Manual inspection of the input file is required.")

parser = argparse.ArgumentParser(description='Correct Newick tree syntax issues')
parser.add_argument('-i', '--input', type=str, required=True, help='Path to input Newick tree file')
parser.add_argument('-o', '--output', type=str, required=True, help='Path to corrected output Newick tree file')
args = parser.parse_args()

correct_newick_tree(args.input, args.output)
