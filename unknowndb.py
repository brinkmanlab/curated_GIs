#!/usr/bin/env python3
"""
Generate unknown gene database for prokka from high identity blast matches
See genes.sh
"""
from Bio import SeqIO, SeqRecord

clusters = []
with open('blastp_genes_filtered.tabular') as f:
    for line in f:
        r = line.split("\t")
        for c in clusters:
            if r[0] in c or r[1] in c:
                c.add(r[0])
                c.add(r[1])
                break
        else:
            clusters.append({r[0], r[1]})

i = 0
while i < len(clusters):
    a = i + 1
    while a < len(clusters):
        if clusters[a] and clusters[i] & clusters[a]:
            clusters[i] |= clusters[a]
            del clusters[a]
        else:
            a += 1
    i += 1

with open('gene_database.faa') as f, open('unknowndb.faa', 'w') as o:
    SeqIO.write((SeqRecord.SeqRecord(id=f'unk{i}.{a}', seq=r.seq, description=f'~~~unk{i}~~~unk{i}~~~') for a, r in enumerate(SeqIO.parse(f, 'fasta')) for i, c in enumerate(clusters) if r.id in c), o, 'fasta')
