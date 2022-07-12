#!/usr/bin/env bash
. ../../deps.sh

biopython.convert -s genomes.gb genbank genome.gb genbank

qlength=$(biopython.convert -q '[length([0].seq)]' query.gb genbank /dev/stdout text)

awk 'BEGIN{OFS="\t"}$5>10 { if ($9 < $10) print $2, $9, $10; else print $2, $10, $9 }' blastn_coverage.tabular | sort -n |
bedtools merge -d 14000 -i /dev/stdin |
awk -v "qlength=$qlength" 'BEGIN{OFS=" "} $3 - $2 > alen[$1] { alen[$1]=$3-$2; best[$1]=($2-(qlength/2) ":" $3+(qlength/2)); } END { for (key in best) print i++, key, best[key]; }' |
xargs -l bash -c "$bpc -q \"[? id == \\\`\\\"\$1\\\"\\\`] | [[0][\$2]]\" genomes.gb genbank \"GI_\${0}.gb\" genbank"

docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/mnt" -w /mnt quay.io/biocontainers/clinker-py:0.0.23--pyh5e36f6f_0 bash -c 'clinker -j 10 -p clinker.html query.gb GI_*.gb'