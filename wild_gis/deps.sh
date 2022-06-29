#!/usr/bin/env bash
shopt -s expand_aliases failglob
set -e

bpc="docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt -w /mnt quay.io/biocontainers/biopython.convert:1.3.3--pyh5e36f6f_0 biopython.convert"
alias bedtools="docker run --rm -i --user $(id -u):$(id -g) -v $PWD:/mnt -w /mnt quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2 bedtools"
alias biopython.convert="$bpc"
alias progressiveMauve="docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt -w /mnt quay.io/biocontainers/mauve:2.4.0.snapshot_2015_02_13--hdfd78af_4 progressiveMauve"
alias blast="docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt -w /mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"