import subprocess

# samples
sra_numbers = [
    "SRR12235513", "SRR12235514", "SRR12235515", "SRR12235516"
    ]

# this will download the .sra files
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip /scratch/jorgepm/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
