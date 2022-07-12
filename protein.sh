#!/usr/bin/env bash
# Generate unknown gene database for prokka
# Depends on output of annotate.sh

. ./deps.sh

shopt -u failglob

jpath=$(cat <<'EOF'
[*].let({genome_id: id}, &(features[?type=='CDS' && qualifiers.translation].{
  id: (qualifiers.protein_id[0] || qualifiers.locus_tag[0]),
  seq: qualifiers.translation[0],
  description: join(' ', [genome_id, to_string(sum([location.start, `1`])), to_string(location.end)])
})) | []
EOF
)

rm -f protein_database.faa
for f in sequences/*.gb; do
  biopython.convert -q "$jpath" "$f" genbank /dev/stdout fasta >> protein_database.faa
done

cat protein_database.faa cog-20.fa > protein_database_cogs.faa

# diamond/proteins :
rm -rf diamond/proteins.*
diamond makedb --in protein_database_cogs.faa -d diamond/proteins

# blastp_proteins.tabular :
diamond blastp --threads 10 --fast -k 2 --evalue 1e-09 --outfmt 6 qseqid sseqid bitscore -d diamond/proteins -q protein_database.faa -o blastp_proteins.tabular

# todo only include COGS that exist within GIs
python3 cogdb.py blastp_proteins.tabular cog-20.cog.csv protein_database_cogs.faa cogs.faa

# diamond/cogs :
rm -rf diamond/cogs.*
diamond makedb --in cogs.faa -d diamond/cogs

diamond blastp --threads 10 --fast -k 1 --evalue 1e-09 --outfmt 6 qtitle stitle -d diamond/cogs -q protein_database.faa -o blastp_gi.tabular

sed 's/\t/ /g' blastp_gi.tabular | sort -bt' ' -k'2,2' -k'3n,4n' |
awk -f create_clusters.gawk > gene_order.txt

