#!/usr/bin/env bash
# Generate unknown gene database for prokka
# Depends on output of annotate.sh

. ./deps.sh

shopt -u failglob

python3 gene_database.py

# blast/genes :
rm -rf blast/genes.*
blast makeblastdb -in gene_database.fna -input_type fasta -dbtype nucl -title GI_genes -hash_index -out blast/genes > makeblastdb_genes.stdout

# blastn_genes.tabular :
blast blastn -num_threads 10 -max_hsps 1 -out blastn_genes.tabular -outfmt '6 qseqid sseqid qcovs qcovhsp length evalue bitscore' -db blast/genes -query gene_database.fna

# Filter alignments removing <95% identity
awk 'BEGIN {FS=OFS="\t"} $4>=95 && $1!=$2 { print }' blastn_genes.tabular > blastn_genes_filtered.tabular

python3 clusterdb.py

# blast/clusters :
rm -rf blast/clusters.*
blast makeblastdb -in clustered_genes.fna -input_type fasta -dbtype nucl -title GI_genes -hash_index -out blast/clusters > makeblastdb_clusters.stdout