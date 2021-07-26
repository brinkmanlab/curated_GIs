# curated_GIs
Curated database of labelled Genomic Islands found in publications

**If you are cloning this repo with the intention of modifying it, you MUST run ./.sqlite.py after cloning**

# Find similar GI entries by sequence
docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 blastn -num_threads 10 -out /mnt/blastn_coverage.tabular -outfmt '6 qseqid sseqid qcovs qcovus' -db /mnt/blast/GIs -query /mnt/sequences.fna
