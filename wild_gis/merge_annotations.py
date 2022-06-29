from Bio import SeqIO, StreamModeError, SeqFeature, Seq
from collections import defaultdict

with open('filtered.tabular') as f, open('genomes.gbk') as s, open('genomes_annotated.gbk', 'w') as out:
    features = defaultdict(list)
    for l in f:
        asc, hsp, bitscore, start, end = l.split(' ')
        start, end = int(start), int(end)
        if start > end:
            ostart = start
            start = end
            end = ostart
        features[asc].append(SeqFeature.SeqFeature(location=SeqFeature.FeatureLocation(start, end), type='mobile_element', id='GI_alignment', qualifiers={"bitscore": bitscore, "hsp": hsp, "note":"GI_alignment"}))
    records = list(SeqIO.parse(s, 'genbank'))
    for record in records:
        record.features.extend(features[record.id])
    SeqIO.write(records, out, 'genbank')

