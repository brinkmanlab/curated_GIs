#!/usr/bin/env bash
# Generate unknown gene database for prokka
# Depends on output of annotate.sh

. ./deps.sh

shopt -u failglob

rm -f protein_database.faa
for f in sequences/*.gb; do
  biopython.convert -q "$(cat proteins.jmespath)" "$f" genbank /dev/stdout fasta >> protein_database.faa
done

cat protein_database.faa cog-20.fa > protein_database_cogs.faa

# diamond/proteins :
rm -rf diamond/proteins.*
diamond makedb --in protein_database_cogs.faa -d diamond/proteins

# blastp_proteins.tab :
diamond blastp --threads 10 --fast -k 2 --evalue 1e-09 --outfmt 6 qseqid sseqid bitscore -d diamond/proteins -q protein_database.faa -o blastp_proteins.tab

python3 cogdb.py blastp_proteins.tab cog-20.cog.csv protein_database_cogs.faa cogs.faa

# diamond/cogs :
rm -rf diamond/cogs.*
diamond makedb --in cogs.faa -d diamond/cogs

diamond blastp --threads 10 --fast -k 1 --evalue 1e-09 --outfmt 6 qtitle stitle -d diamond/cogs -q protein_database.faa -o blastp_gi.tab

sed 's/\t/ /g' blastp_gi.tab | sort -bt' ' -k'2,2' -k'3n,4n' |
awk -f create_clusters.gawk > gene_order.txt

