#!/usr/bin/env python3
"""
Generate unknown gene database for prokka from high identity blast matches
See genes.sh
"""
import glob

from Bio import SeqIO, SeqRecord

clusters = []
with open('blastn_genes_filtered.tabular') as f:
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

genes = {g: i for i, c in enumerate(clusters) for g in c}

with open('gene_clusters.txt', 'w') as f:
    f.writelines(' '.join(cluster) + '\n' for cluster in clusters)

with open('gene_order.txt', 'w') as order:
    for p in glob.glob('sequences/*.gb'):
        with open(p) as f:
            pairs = []
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
                    if gene_id in genes:
                        gene_id = f"cluster{genes[gene_id]}"
                    if last:
                        pairs.append(f"{last}_{gene_id}" if last < gene_id else f"{gene_id}_{last}")
                    last = gene_id
            if len(pairs) == 0:
                print(p, "has one or no genes")
            else:
                for pair in pairs:
                    print(r.id, pair, file=order)


with open('gene_database.fna') as f, open('clustered_genes.fna', 'w') as o:
    SeqIO.write((SeqRecord.SeqRecord(id=f'cluster{i}', seq=r.seq, description=r.description) for a, r in enumerate(SeqIO.parse(f, 'fasta')) for i, c in enumerate(clusters) if r.id in c), o, 'fasta')
# .{a}
