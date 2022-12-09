#!/bin/bash

# Script by Jorge L. Perez-Moreno, Ph.D.  jorgepm@colostate.edu

fraction="0.1" 	#fraction of reads to keep

# Execute script subsample.sh from base folder containing all treatment reads in corresponding subfolders

for treatment in ./*/ ; do

	cd $treatment
	echo $treatment
	gunzip *.gz
	for file in *_R1_001.fastq ; do
		sample=$(echo $file| cut -d'_' -f 1,3)
		samplein=$(echo $file| cut -d'_' -f 1,2,3)
		python /home/jorgepm/bin/subsample.py $fraction ${samplein}_R1_001.fastq ${samplein}_R2_001.fastq ${sample}_${fraction}_1.fastq ${sample}_${fraction}_2.fastq
#		rm  ${sample}_R1_001.fastq ${sample}_R2_001.fastq
	done
	gzip *.fastq
	cd ..
done


	
