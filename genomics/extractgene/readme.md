blast_pipeline.py is deprecated.

please use extractgene.py instead.

pipeline_utils.py should be in the same folder as extractgene.py.

extractgene.py will 1) index your genome using hisat2 OR create a BLAST db, 2) map your transcript OR blast your DNA sequence, 3) extract a sequence of X length upstream of the first exon / first BLAST HSP
