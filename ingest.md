# Files 
- blast/GIs : `docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 makeblastdb -in /mnt/sequences.fna -input_type fasta -dbtype nucl -title GIs -hash_index -out /mnt/blast/GIs > makeblastdb.stdout`
- sequences.fna : `awk 1 sequences/* > sequences.fna`
- blastn_coverage.tabular : `docker run --rm --user $(id -u):$(id -g) -v $PWD:/mnt quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 blastn -num_threads 10 -out /mnt/blastn_coverage.tabular -outfmt '6 qseqid sseqid qcovs qcovus' -db /mnt/blast/GIs -query /mnt/sequences.fna`
- blastn_coverage_filtered.tabular : `awk 'BEGIN {FS=OFS="\t"} $4>=30 && $1!=$2 { print }' blastn_coverage.tabular > blastn_coverage_filtered.tabular`



# GI Curation Troubleshooting

## Alternate Accessions for a Strain
Problem: Some strains have been sequenced multiple times, and therefore have multiple valid NCBI accessions associated with the same strain name

Solution: No longer deduplicate in cases where the same strain name has alternate accessions values associated with it - added as distinct entries in the Strain table if the accessions are unique.

## Strain without Accession
Problem: Some strains do not have any accession value associated with them (this is often the case when the sequence for the GI/phage has its own accession and referencing the sequence of the strain is rendered unnecessary)

Solution: The database was altered slightly to allow for a null value associated with a strain. (Given the problem discussed above related to alternate accessions for a strain, it is possible to have multiple entries for strains with e same name, some with an accession/gbuid and one with a null value)

## Ambiguous GI Names

Problem: Some GIs are given uninformative names, often GI1,...,GI<n> where GI<n> is the nth GI discovered and characterized for that species. This can lead to ambiguous and confusing cases (eg. GI1 from Enterococcus faecium is not the same as GI1 from Burkholderia cenocepacia) and if these names are not distinguished in the output the user could reach inaccurate conclusions

Solution: Duplicate GIs in the database will have a name that combines their given name in the literature and a unique curated genomic island (cGI) ID assigned at the time of database development (ex. GI1_cGI_42).

## Multiple Publications per Entry

Problem: Multiple publications were appended and separated by semicolons (usually) in the input data. These were being stored in a single field.

Solution: A separate publications table was generated to store the publication information so that each publication can have a unique entry. This is in keeping with standard database development practices.

## Characterized GIs without Names

Problem: Some GIs that were well studied and characterized were not assigned names by the publishing authors, but were still collected in the curation process.

Solution: These GIs will only have the ‘cGI#’ ID for their name

## Multiple Chromosomes Associated with a Genome/Strain

Problem: Each chromosome from a genome with multiple chromosomes will have a unique NCBI accession/entry associated with the same strain name. These were originally being deduplicated.

Solution: This problem will be solved by the implementation of the solution for the “Alternate Accessions for a Strain” solution, which now allows entries with the same strain name to be entered separately and keep their alternate accession IDs

## Accessions for phage/GIs versus genomes

Problem: In some cases, the accession associated with an entry points to a genome (from which the GI sequence spans a given set of coordinates) and in other cases the accession value points to an entry only containing the phage/GI sequence itself

Solution: In the case where the accession is associated with just the phage/GI sequence, the Sequence entry for that GI will have a gbuid that points directly to the NCBI entry (as opposed to the strain).

## Missing sequences

Problem: A record is missing a sequence fasta file

Solution: The GI accession or strain accession is searched and downloaded from the NCBI Nucleotide database. If the strain accession is used, the start and end coordinates must be present and only that coordinate range is retained in the fasta

## GI deduplication

Problem: Some records refer to the same GI

Solution: If records have the same GI accession, or strain accession with same start/end coordinates, then the GI will be deduplicated.