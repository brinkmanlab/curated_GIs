#!/usr/bin/env python3
"""
Generate predicted gene database for blast to do an all vs all alignment
See genes.sh
"""

from Bio import SeqIO, SeqRecord
import glob
from collections import defaultdict

genes = defaultdict(int)

with open('gene_database.fna', 'w') as out: #, open('gene_order.txt', 'w') as order:
    for fn in glob.glob('sequences/*.gb'):
        #pairs = []
        print(fn)
        with open(fn) as f:
            r = next(SeqIO.parse(f, 'genbank'))
            last = None
            for feat in r.features:
                if feat.type == 'CDS':
                    gene_id = feat.qualifiers.get('gene')
                    if not gene_id:
                        gene_id = feat.qualifiers.get('note')
                        if gene_id:
                            gene_id = gene_id[0].split(';')
                            if len(gene_id) > 1 and ' ' not in gene_id[0]:
                                gene_id = gene_id[0]
                            else:
                                gene_id = None
                    if not gene_id:
                        gene_id = feat.qualifiers.get('locus_tag')
                    if not gene_id:
                        gene_id = feat.qualifiers.get('protein_id')
                    if isinstance(gene_id, list): gene_id = gene_id[0]
                    gene_id = gene_id.replace(' ', '_')
                    #if last:
                    #    pairs.append(f"{last}_{gene_id}" if last < gene_id else f"{gene_id}_{last}")
                    #last = gene_id
                    #genes[gene_id] += 1
                    #if genes[gene_id] > 1:
                    #    gene_id += f"_{genes[gene_id]}"
                    SeqIO.write(SeqRecord.SeqRecord(seq=feat.extract(r.seq), id=gene_id, description=fn), out, 'fasta')
            #print(r.id, ' '.join(pairs), file=order)
