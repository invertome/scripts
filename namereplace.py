# replace dammit! transcript names for original Trinity contig IDs using namemap
# script by Jorge L. Perez-Moreno, Ph.D.
# usage: python namereplace.py namemap_file input output

import sys

dictionary = sys.argv[1]

origin = sys.argv[2]

results = sys.argv[3]

# opening csv file, read mode, comma as delimeter

import csv
reader = csv.reader(open(dictionary, 'rb'), delimiter=',')

# creating a dictionary using each row as a combination of key, value

mydict = dict((rows[1],rows[0]) for rows in reader)
#print(mydict)

# writing newly formed dictionary into a new file, preserving it's object type with pickle
import pickle

with open('mapping.csv', 'wb') as f:
	pickle.dump(mydict, f)

# import dependencies

import re

# open dictionary as python object
with open('mapping.csv', 'rb') as newdict:
	mydict = pickle.loads(newdict.read())


# define replacing functions
def replacer_factory(spelling_dict):
	def replacer(match):
		word = match.group()
		return spelling_dict.get(word, word)
	return replacer

def reps(text):
	pattern = r'\w+(?=\.p)'  # pattern to match exact IDs and not similar ones
	replacer = replacer_factory(mydict)
	return re.sub(pattern, replacer, text)

def main():
# open original peptide file as in file
	with open(origin) as in_file:
		text = in_file.read()
# open new results peptide file as out file
	with open(results, 'w') as out_file:
		out_file.write(reps(text))
# print results on screen to verify
#	print(reps(text))


if __name__ == '__main__':
	main()


print('DONE! Results saved as: ' + results)

