#!/usr/bin/env bash

. ../../deps.sh

prokka --force --compliant --proteins ../../cogs.faa --outdir prokka/ genomes.fna

#  makeblastdb -dbtype prot -in \.\.\/\.\.\/cogs\.faa -out prokka\/\/proteins -logfile /dev/null

# cat prokka\/\/PROKKA_07122022\.proteins\.tmp\.1\.faa |
# parallel --gnu --plain -j 8 --block 922001 --recstart '>' --pipe blastp -query - -db prokka//proteins -evalue 1e-09 -qcov_hsp_perc 80 -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no > prokka\/\/PROKKA_07122022\.proteins\.tmp\.1\.blast 2> /dev/null