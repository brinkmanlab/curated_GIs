#!/usr/bin/env bash
. ../../deps.sh

cutoff=0.5

biopython.convert -q "$(cat ../../proteins.jmespath)" genomes.gb genbank genomes.faa fasta

biopython.convert -q "$(cat ../../CDS.jmespath)" genomes.gb genbank cds.tab text

#blast blastp -evalue 1e-09 -qcov_hsp_perc 95 -num_threads 10 -max_hsps 1 -num_alignments 1 -seg no -outfmt '6 qseqid sseqid qcovs qcovus qcovhsp length evalue bitscore' -db ../../blast/cogs -query genomes.faa -out blast.tab
diamond blastp --threads 10 --fast -k 1 --evalue 1e-09 --outfmt 6 qtitle stitle -d ../../diamond/cogs -q genomes.faa -o blast.tab

sed 's/\t/ /g' blast.tab cds.tab | sort -bt' ' -k'2,2' -k'3n,4n' > cds_combined.tab
awk -v cutoff=$cutoff -f ../create_set.gawk cds_combined.tab | sort | uniq > blast_sets.tab

setsimilaritysearch --input-sets ../../gene_order.txt blast_sets.tab --output-pairs jaccard.csv --similarity-threshold cutoff
awk 'BEGIN{ FS=","; OFS="\t"}NR>1{ match($1, "([^:]+):([^-]+)-(.+)", a); print a[1], a[2], a[3], $2, $5 }' jaccard.csv > jaccard.bed


# Output summary
matches=$(cut -f1 jaccard.bed | uniq | wc -l)

biopython.convert -q '[*].id' genomes.gb genbank genomes.ids text
genomes=$(wc -l < genomes.ids)

echo "$matches genomes with GI detected of ${genomes}. No GIs detected in:" > report.txt

cat <(cut -f1 -d$'\t' jaccard.bed | uniq) genomes.ids | sort | uniq -u >> report.txt

cat report.txt

# awk '$1 == $4' jaccard.bed > recall.bed
# awk '$1 != $4 && $5 == 1' jaccard.bed > identical.bed