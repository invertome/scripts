# script to subsample large paired-end fasta files (i.e. Illumina reads)
# called by: python subsample.py <fraction> <input file 1> <input file 2> <output file 1> <output file 2>
# fraction is a percentage, from 0 to 1


import sys, random, itertools
import HTSeq

fraction = float( sys.argv[1] )
in1 = iter( HTSeq.FastqReader( sys.argv[2] ) )
in2 = iter( HTSeq.FastqReader( sys.argv[3] ) )
out1 = open( sys.argv[4], "w" )
out2 = open( sys.argv[5], "w" )

for read1, read2 in zip( in1, in2 ):
   if random.random() < fraction:
      read1.write_to_fastq_file( out1 )
      read2.write_to_fastq_file( out2 )
      
out1.close()
out2.close()
