#!/usr/bin/env python3
"""
Generate database for prokka from high identity blast matches
See protein.sh
"""

import sys
import csv

from Bio import SeqIO

_, blast_output_path, cog20_index_path, protein_database_cogs, output_path = sys.argv

# Columns:
# 1.	Gene ID (GenBank or ad hoc)
# 2.	NCBI Assembly ID
# 3.	Protein ID (GenBank if conforms to [A-Za-z0-9_]+\.[0-9]+ regex; ad hoc otherwise)
# 4.	Protein length
# 5.	COG footprint coordinates on the protein. "201-400" means "from position 201 to position 400"; "1-100=201-300" indicates a segmented footprint, 1-100 AND 201-300
# 6.	Length of the COG footprint on the proteins
# 7.	COG ID
# 8.	reserved
# 9.	COG membership class (0: footprint covers most of the protein and most of the COG profile; 1: footprint covers most of the COG profile and part of the protein; 2: footprint covers most of the protein and part of the COG profile; 3: partial match on both protein and COG profile)
# 10.	PSI-BLAST bit score for the match between the protein and COG profile
# 11.	PSI-BLAST e-value for the match between the protein and COG profile
# 12.	COG profile length
# 13.	Protein footprint coordinates on the COG profile
with open(cog20_index_path) as f:
    cog20_index = {pid.replace('.', '_'): (cogid, gene) for gene, _, pid, _, _, _, cogid, *_ in csv.reader(f)}

# Find anything that aligns to a cog20 protein and assign it the same COG id
cog_hits = {k: (sys.maxsize, v) for k, v in cog20_index.items()}
non_terminal_hits = {}
with open(blast_output_path) as f:
    for dbid, qid, bitscore in csv.reader(f, 'excel-tab'):
        qid = qid.replace('.', '_')
        dbid = dbid.replace('.', '_')
        bitscore = float(bitscore)
        if dbid == qid:
            continue
        hit = cog_hits.get(qid)
        if not hit or hit[0] < bitscore:
            cog = cog20_index.get(dbid)
            if cog:
                cog_hits[qid] = (bitscore, cog[0])
            else:
                hit = non_terminal_hits.get(qid)
                if not hit or hit[0] < bitscore:
                    non_terminal_hits[qid] = (bitscore, dbid)

# Cluster anything that did not point directly to a COG protein
clusters = []
for qid, (_, dbid) in non_terminal_hits.items():
    for c in clusters:
        if qid in c or dbid in c:
            c.add(qid)
            c.add(dbid)
            break
    else:
        clusters.append({qid, dbid})

print("Collapsing clusters")
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

print("Clustered")

# Use any second order alignments to a COG protein to assign the entire cluster to the same COG
non_terminal_hits.update({v: (b, k) for k, (b, v) in non_terminal_hits.items()})
keys = set(cog_hits.keys())
todel = []
for i, c in enumerate(clusters):
    cluster_hits = c.intersection(keys)
    if len(cluster_hits) == 0: continue
    chit = next(h for h in cluster_hits)
    if len(cluster_hits) > 1:
        hit = non_terminal_hits[chit]
        for ch in cluster_hits:
            newhit = non_terminal_hits[ch]
            if hit[0] < newhit[0]:
                hit = newhit
                chit = ch
    for p in c:
        cog_hits[p] = cog_hits[chit]
    todel.append(i)

print("Clusters associated")

for d in reversed(todel):
    del clusters[d]

for i, c in enumerate(clusters):
    for p in c:
        cog_hits[p] = (None, f"pCOG{i}")

cog_hits = {k: v[1] for k, v in cog_hits.items()}
genes = {k: v for (k, v) in cog20_index.values()}

with open(protein_database_cogs) as f, open(output_path, 'w') as o:
    for r in SeqIO.parse(f, 'fasta'):
        r.description = r.description.lstrip(r.id + " ")
        r.id = r.id.replace('.', '_')
        gene = ''
        desc = ''
        if '~~~' not in r.description:
            desc = r.description.rsplit('[', 1)[0].strip(' ')
        cogid = ''
        cog = cog20_index.get(r.id)
        if cog:
            cogid, gene = cog
        else:
            cogid = cog_hits.get(r.id, '')
        if not gene:
            gene = genes.get(cogid, '')
        if '~~~' in r.description:
            descs = r.description.split('~~~')
            if not gene and descs[1] != 'None':
                gene = descs[1]
            if not desc and descs[2] != 'None':
                desc = descs[2]

        #r.description = f"~~~{gene}~~~{desc}~~~{cogid}"  # for prokka db
        r.description = cogid or f"{r.id} {r.id}"  # SeqIO will trim first copy
        SeqIO.write(r, o, 'fasta')