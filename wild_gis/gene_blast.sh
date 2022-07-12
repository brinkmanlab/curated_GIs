#!/usr/bin/env bash
. ../../deps.sh

jpath=$(cat <<'EOF'
[*].let({genome_id: id}, &(features[?type=='CDS' && qualifiers.translation].{
  id: join(' ', [genome_id, (qualifiers.protein_id[0] || qualifiers.locus_tag[0])]),
  seq: qualifiers.translation[0],
  description: join('', [to_string(sum([location.start, `1`])), ' ', to_string(location.end)])
})) | []
EOF
)

biopython.convert -q "$jpath" genomes.gb genbank genomes.faa fasta

#blast blastp -evalue 1e-09 -qcov_hsp_perc 95 -num_threads 10 -max_hsps 1 -num_alignments 1 -seg no -outfmt '6 qseqid sseqid qcovs qcovus qcovhsp length evalue bitscore' -db ../../blast/cogs -query genomes.faa -out genes_blast.tabular
diamond blastp --threads 10 --fast -k 1 --evalue 1e-09 --outfmt 6 qtitle stitle -d ../../diamond/cogs -q genomes.faa -o genes_blast.tabular

sed 's/\t/ /g' genes_blast.tabular | sort -bt' ' -k'1,1' -k'3n,4n' |
awk -f ../create_set.gawk > genes_blast_sets.tabular

setsimilaritysearch --input-sets ../../gene_order.txt genes_blast_sets.tabular --output-pairs jaccard.csv --similarity-func jaccard --similarity-threshold 0.5
awk 'BEGIN{ FS=","; OFS="\t"}NR>1{ match($1, "([^:]+):(.+)", a); print a[1], a[2], a[2]+1000, $2, $5 }' jaccard.csv > jaccard.bed

# intersect=$(wc -l $(uniq -d pairs))
# union=$(wc -l $(uniq pairs))
# Jacard=(( $intersect / $union ))