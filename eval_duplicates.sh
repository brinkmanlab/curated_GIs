#!/usr/bin/env bash
# Align all vs all sequences generating blastn_alignments.html
set -ev

# sequences.fna :
awk 1 sequences/*.fasta sequences/*.fna > sequences.fna

sed -i 's/φ/phi/' sequences.fna

# blast/GIs :
rm -rf blast/GIs.*
docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 makeblastdb -in /mnt/sequences.fna -input_type fasta -dbtype nucl -title GIs -hash_index -out /mnt/blast/GIs > makeblastdb.stdout

# blastn_coverage.tabular : 
docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 blastn -num_threads 10 -out /mnt/blastn_coverage.tabular -outfmt '6 qseqid sseqid qcovs qcovus qcovhsp length evalue bitscore' -db /mnt/blast/GIs -query /mnt/sequences.fna

# blastn_coverage_filtered.tabular : 
awk 'BEGIN {FS=OFS="\t"} $4>=30 && $1!=$2 { print }' blastn_coverage.tabular > blastn_coverage_filtered.tabular

# blastn_coverage_filtered_paths.tabular : Import blastn_coverage_filtered.tabular into temp table and execute 
sqlite3 GIs.sqlite <<EOF > blastn_coverage_filtered_paths.tabular
CREATE TEMP TABLE blastn_coverage_filtered ( C1 text, C2 text, C3 text, C4 text, qcovhsp integer, length integer, evalue text, bitscore integer );
.mode tabs
.import blastn_coverage_filtered.tabular blastn_coverage_filtered
SELECT seq1.id, seq2.id, f.C3, f.C4, seq1.path, seq2.path, f.qcovhsp, f.length, f.evalue, f.bitscore from blastn_coverage_filtered as f, sequences as seq1, sequences as seq2 where f.C1 like '%|' || seq1.id and f.C2 like '%|' || seq2.id and seq1.gi != seq2.gi;
EOF

# blastn_coverage_filtered_paths_dedup.tabular :
# TODO dont dedup, reverse pair should have same coverage unless one is inside the other. Keep reverse qcov in report
awk '{if ($1 > $2) {key = $1$2;} else {key = $2$1;}; if (dup[key] != 1) {print;}; dup[key] = 1;}' blastn_coverage_filtered_paths.tabular > blastn_coverage_filtered_paths_dedup.tabular

rm -rf alignments/*
# alignments/*.xmfa : 
cut -f 1,2,5,6 blastn_coverage_filtered_paths_dedup.tabular | while IFS=$'\t' read l r lp rp; do
  docker run --rm --user $(id -u):$(id -g) -v "$PWD:/mnt" quay.io/biocontainers/mauve:2.4.0.snapshot_2015_02_13--hdfd78af_4 progressiveMauve --output=/mnt/alignments/$l-$r.xfma /mnt/$lp /mnt/$rp
done

# alignments/*.png : 
# Mauve crashes half way through, likely a memory leak. This step needs to just be restarted, the macro will resume for only the missing pictures.
cid=$(docker run -d --rm --user $(id -u):$(id -g) -v "$HOME/.Xauthority:$HOME/.Xauthority:rw" --env="XAUTHORITY=$HOME/.Xauthority" --env="DISPLAY" -v "$PWD:/mnt" --net=host quay.io/biocontainers/mauve:2.4.0.snapshot_2015_02_13--hdfd78af_4 Mauve)
sleep 2
actiona -stex ./gen_images.ascr
docker kill $cid

# blastn_coverage_filtered_paths_dedup_pic.tabular : 
sqlite3 GIs.sqlite <<EOF > blastn_coverage_filtered_paths_dedup_pic.tabular
CREATE TEMP TABLE blastn_coverage_filtered ( seq1 text, seq2 text, C3 text, C4 text, path1 text, path2 text );
.mode tabs
.import blastn_coverage_filtered_paths_dedup.tabular blastn_coverage_filtered
SELECT gi1.name, gi1.id, seq1.id, str1.name, gi2.name, gi2.id, seq2.id, str2.name, f.C3, f.C4, 'alignments/' || seq1.id || '-' || seq2.id || '.xfma.png' as pic
from blastn_coverage_filtered as f, sources as src1, sources as src2, strains as str1, strains as str2, sequences as seq1, sequences as seq2, genomic_islands as gi1, genomic_islands as gi2
where seq1.id = f.seq1 and seq2.id = f.seq2 and str1.id = src1.strain and seq1.id = src1.seq and str2.id = src2.strain and seq2.id = src2.seq and gi1.id = seq1.gi and gi2.id = seq2.gi;
EOF

cat <<EOF > blastn_aligments.html
<html><body>
<style>td { border-bottom: 1px solid black;} img { height: 397px; }</style>
<table>
<tr><td>C1<br/>name</td><td>C2<br/>name</td><td>qcovs</td><td>qcovus</td><td></td></tr>
EOF

# blastn_alignments.html : 
awk 'BEGIN {FS="\t";}{ print "<tr><td>"$1"|"$2"|"$3"<br/>"$4"</td><td>"$5"|"$6"|"$7"<br/>"$8"</td><td>"$9"</td><td>"$10"</td><td><img src=\""$11"\" /></td></tr>\n"}' blastn_coverage_filtered_paths_dedup_pic.tabular >> blastn_aligments.html

cat <<EOF >> blastn_aligments.html
</table>
</body></html>
EOF