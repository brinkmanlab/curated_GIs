#!/usr/bin/env bash

shopt -s expand_aliases failglob
set -e

export PROJECT=$(dirname $(realpath "$BASH_SOURCE"))
wd=$(realpath --relative-to="$PROJECT" "$PWD")

bpc="docker run --rm --user $(id -u):$(id -g) -v $PROJECT:/mnt -w /mnt/$wd quay.io/biocontainers/biopython.convert:1.3.3--pyh5e36f6f_0 biopython.convert"
alias bedtools="docker run --rm -i --user $(id -u):$(id -g) -v $PROJECT:/mnt -w /mnt/$wd quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2 bedtools"
alias biopython.convert="$bpc"
alias progressiveMauve="docker run --rm --user $(id -u):$(id -g) -v $PROJECT:/mnt -w /mnt/$wd quay.io/biocontainers/mauve:2.4.0.snapshot_2015_02_13--hdfd78af_4 progressiveMauve"
alias blast="docker run --rm --user $(id -u):$(id -g) -v $PROJECT:/mnt -w /mnt/$wd quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
alias setsimilaritysearch="docker run --rm --user $(id -u):$(id -g) -v $PROJECT:/mnt -w /mnt/$wd quay.io/biocontainers/setsimilaritysearch:1.0.0 all_pairs.py"
alias prokka="docker run --rm --user $(id -u):$(id -g) -v$PROJECT:/mnt -w /mnt/$wd staphb/prokka:latest prokka"
alias diamond="docker run --rm --user $(id -u):$(id -g) -v$PROJECT:/mnt -w /mnt/$wd quay.io/biocontainers/diamond:2.0.15--hb97b32f_0 diamond"