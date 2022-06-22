#!/usr/bin/env bash
# Generate unknown gene database for prokka
# Depends on output of annotate.sh

python3 gene_database.py

# blast/genes :
rm -rf blast/genes.*
docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 makeblastdb -in /mnt/gene_database.faa -input_type fasta -dbtype prot -title GI_genes -hash_index -out /mnt/blast/genes > makeblastdb_genes.stdout

# blastp_genes.tabular :
docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 blastp -num_threads 10 -max_hsps 1 -out /mnt/blastp_genes.tabular -outfmt '6 qseqid sseqid qcovs qcovhsp length evalue bitscore' -db /mnt/blast/genes -query /mnt/gene_database.faa

# Filter alignments removing <95% identity
awk 'BEGIN {FS=OFS="\t"} $4>=95 && $1!=$2 { print }' blastp_genes.tabular > blastp_genes_filtered.tabular

python3 unknowndb.py