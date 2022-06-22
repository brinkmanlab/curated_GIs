#!/usr/bin/env bash
# Annotate all sequences with predicted genes
#shopt -s failglob

mkdir -p genes
rm -f genes/*
#for f in sequences/*.fna; do
#  n=$(basename "$f" '.fna')
#  # https://github.com/tseemann/prokka#command-line-options
#  docker run --rm --user $(id -u):$(id -g) -v$PWD:/mnt staphb/prokka:latest prokka --force --compliant --proteins /mnt/unknowndb.faa --outdir /mnt/prokka/ --prefix $n /mnt/$f
#done

for f in sequences/*.gb; do
  docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/biopython.convert:1.3.1--pyh5e36f6f_0 biopython.convert /mnt/${f} genbank /mnt/sequences/$(basename $f .gb).gff gff
done

sed -i '1 i\##gff-version 3' sequences/*.gff

for f in sequences/*.gff; do
  docker run --rm --user $(id -u):$(id -g) -v/usr/share/fonts:/usr/share/fonts -v$PWD:/mnt quay.io/biocontainers/genometools-genometools:1.6.1--py38h23571c4_2 gt sketch -v -flattenfiles no --width 1400 --format svg "/mnt/genes/$(basename $f '.gff').svg" /mnt/$f
done