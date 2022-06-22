#!/usr/bin/env python3
"""
Generate predicted gene database for blast to do an all vs all alignment
See genes.sh
"""

from Bio import SeqIO, SeqRecord, Seq
import glob

with open('gene_database.faa', 'w') as out:
    for fn in glob.glob('prokka/*.gbk'):
        print(fn)
        with open(fn) as f:
            for r in SeqIO.parse(f, 'genbank'):
                SeqIO.write((SeqRecord.SeqRecord(seq=Seq.Seq(feat.qualifiers['translation'][0]), id=feat.qualifiers['locus_tag'][0], description=fn) for feat in r.features if feat.type == 'CDS' and 'gene' not in feat.qualifiers), out, 'fasta')
