#!/usr/bin/env python3

import sys
from Bio import SeqIO
from mimetypes import guess_type
import gzip

filename = sys.argv[1]

encoding = guess_type(filename)[1]
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

with _open(filename) as ffile:
    for seq_record in SeqIO.parse(ffile, 'fasta'):
        if seq_record.seq.count('X') < 1:
            print(seq_record.format("fasta"))
