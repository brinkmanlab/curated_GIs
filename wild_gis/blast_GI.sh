#!/usr/bin/env bash
. ../deps.sh

biopython.convert genomes.gb genbank genomes.fna fasta
biopython.convert query.gb genbank query.fna fasta

rm -rf blast/GIs.*
blast makeblastdb -in genomes.fna -input_type fasta -dbtype nucl -title GIs -hash_index -out blast/GIs > makeblastdb.stdout

blast blastn -num_threads 10 -out blastn_coverage.tabular -outfmt '6 qseqid sseqid qcovs qcovus qcovhsp length evalue bitscore sstart send' -db blast/GIs -query query.fna

awk '$5>10 { print $2, $5, $6, $9, $10 }' blastn_coverage.tabular > filtered.tabular

#awk -F' ' '{ OFS="\t"; print $1, "blast", "genomic_island", $4, $5, $2, ".",".","bitscore="$3}' filtered.tabular > alignments.gff3

#python ../merge_annotations.py