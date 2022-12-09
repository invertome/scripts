# Import necessary modules
import argparse
from ete3 import Tree
import statistics

# Parse the command-line arguments
parser = argparse.ArgumentParser(description='Prune long branches from a phylogenetic tree in Newick format.')
parser.add_argument('-in', '--input_file', required=True, help='the input Newick format file')
parser.add_argument('-out', '--output_file', required=True, help='the output Newick format file')
parser.add_argument('-pruned', '--summary_file', required=True, help='the output file for the summary of pruned taxa')
parser.add_argument('-n', '--num_stddev', type=int, default=3, help='the number of standard deviations used to determine whether a branch is long')
args = parser.parse_args()

# Load the phylogenetic tree from a Newick format file
tree = Tree(args.input_file)

# Define a function to calculate the mean branch length within a clade
def mean_branch_length(clade):
  # Initialize a list to store the branch lengths
  branch_lengths = []

  # Iterate over the clade's descendents
  for descendent in clade.get_descendants():
    # Add the descendent's branch length to the list
    branch_lengths.append(descendent.dist)

  # Return the mean branch length for the clade
  return statistics.mean(branch_lengths)

# Define a function to prune long branches from the phylogenetic tree
def prune_long_branches(tree, num_stddev):
  # Iterate over the clades in the phylogenetic tree
  for clade in tree.get_descendants():
    # Calculate the mean branch length for the clade
    mean_length = mean_branch_length(clade)

    # Calculate the standard deviation of the branch lengths for the clade
    stddev = statistics.stdev(clade.get_descendants(), mean_length)

    # Iterate over the clade's descendents
    for descendent in clade.get_descendants():
      # If the descendent's branch length exceeds 3 times the standard deviation of the mean branch length within the clade, prune the descendent from the tree
      if descendent.dist > num_stddev * stddev:
        descendent.delete()

# Prune long branches from the phylogenetic tree
prune_long_branches(tree)

# Output the pruned phylogenetic tree to a Newick format file
tree.write(format=1, outfile=args.output_file)

# Output a summary of the pruned taxa to a text file
with open(args.summary_file, 'w') as f:
  for node in tree.traverse():
    if node.is_leaf() and node.is_deleted():
      f.write(node.name + '\n')
