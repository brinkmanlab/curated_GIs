# curated_GIs
Curated database of labelled Genomic Islands found in publications

**If you are cloning this repo with the intention of modifying it, you MUST run ./.sqlite.py after cloning**


BLAST Analysis:
qcovs - the percent of no. of bases in the query sequence aligned with the subject sequence (match or mismatch). 
    The bases can be in one HSP or several HSPs (overlap) but they are counted only once. Gaps (in the subject sequence) are treated as mismatches.
qcovus - is a measure of Query Coverage that counts a position in a subject sequence for this measure only once. 
    The second time the position is aligned to the query is not counted towards this measure.


GI genes clustered, including:
https://www.ncbi.nlm.nih.gov/research/cog/