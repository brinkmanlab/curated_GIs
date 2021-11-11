import csv
import sqlite3
import re
import os
import sys
import pickle
from collections import defaultdict
from urllib.error import HTTPError

from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC

Entrez.email = 'nolan_w@sfu.ca'

stdout = sys.stdout
sys.stdout = open('ingest.log', 'w')
os.remove('GIs.sqlite')
conn = sqlite3.connect('GIs.sqlite')
conn.row_factory = sqlite3.Row
with open('schema.sql') as f:
    conn.executescript(f.read())

strain_names = {}
if os.path.isfile('strain_names.dump'):
    with open('strain_names.dump', 'rb') as f:
        strain_names = pickle.load(f)


def fetch(line, gbuid, path, slice=None, id=None, description=None):
    try:
        handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=gbuid)
        records = list(SeqIO.parse(handle, 'fasta'))
        if len(records) < 1:
            print(line, gbuid, "No sequence returned by NCBI")
            return None
        if len(records) > 1:
            print(line, gbuid, "More than one sequence returned by NCBI", tuple(r.id for r in records))
            return None
        if len(records[0].seq) == 0:
            print(line, gbuid, "Sequence returned by NCBI is empty", records[0])
            return None
        if slice is not None:
            if slice.start >= slice.stop:
                # TODO compliment sequence?
                records[0] = records[0][slice.stop:slice.start]
            else:
                records[0] = records[0][slice]
        if id is not None:
            records[0].id = id
        if description is not None:
            records[0].description = description
        # scrub bad chars from description
        records[0].description = records[0].description.replace('ﬄ', 'ffl').replace('\u2010', '-')
        with open(path, 'w') as file:
            SeqIO.write(records, file, 'fasta')
        return round(GC(records[0].seq), 1)
    except HTTPError:
        print(line, gbuid, "No sequence returned by NCBI or bad request")
        return None


def strain_name(line, gbuid, name):
    if not gbuid:
        return name
    new_name = ''
    if name in strain_names:
        new_name = strain_names[gbuid]
    else:
        try:
            handle = Entrez.esummary(db="nucleotide", id=gbuid)
            records = Entrez.read(handle)
            handle.close()
            if len(records) < 1:
                print(line, gbuid, "Strain gbuid returned no results by NCBI")
                return name
            if len(records) > 1:
                print(line, gbuid, "More than one record returned by NCBI", tuple(r.id for r in records))
                return name
            taxid = records[0]['TaxId']
            handle = Entrez.efetch(db="taxonomy", id=taxid)
            records = Entrez.read(handle)
            handle.close()
            if len(records) < 1:
                print(line, gbuid, "Strain taxid returned no results by NCBI", taxid)
                return name
            if len(records) > 1:
                print(line, gbuid, "More than one record returned by NCBI for taxid", taxid, tuple(r.id for r in records))
                return name
            new_name = records[0]['ScientificName']
            strain_names[gbuid] = new_name
        except HTTPError:
            print(line, gbuid, "No records returned by NCBI or bad request")
            return name
    if name != new_name:
        print(line, gbuid, f"Provided strain name ({name}) replaced with NCBI strain name ({records[0]['ScientificName']})")
    return new_name


def dup(line, table, new, old, quiet=False):
    if old is None:
        if not quiet:
            print(f'{line} {table}: failed to meet constraint, discarding {tuple(new.values())}')
        return None
    elif any(old[i] != new[i] for i in new.keys()):
        print(f'{line} {table}: duplicate  {tuple(old or (old,))}')
        print(f'{line} {table}: discarding {tuple(new.values())}')
    return old['id']


with open('GIs.csv') as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    line = 1
    gis = defaultdict(int)
    for row in reader:
        line += 1
        # genomic_islands
        gi = None
        name = row[0].strip().replace(' ', ' ')
        if not name:
            print(line, "No GI name, substituting 'cGI'", row)
            name = "cGI"
        type = row[1].strip().replace(' ', ' ') or None
        role = row[2].strip().replace(' ', ' ') or None
        # TODO deduplicate these variables from below
        strain_gbuid = row[5].strip() or None # TODO handle multiple accessions
        gi_gbuid = row[6].strip() or None
        try: start = int(row[8]) if row[8] else None
        except: start = None
        try: end = int(row[9]) if row[9] else None
        except: end = None
        vals = dict(name=name, type=type, role=role)
        # deduplicate gi on source accession and start-end or gi accession
        if gi_gbuid:
            gi = dup(line, 'genomic_islands', vals, conn.execute(f'''SELECT genomic_islands.* FROM genomic_islands, sequences WHERE gbuid = ? AND genomic_islands.id = sequences.gi''', (gi_gbuid,)).fetchone(), True)
        if not gi and strain_gbuid and start is not None and end is not None:
            gi = dup(line, 'genomic_islands', vals, conn.execute(f'''
                SELECT g.* FROM genomic_islands as g, sources as src, sequences as seq, strains as str 
                WHERE src.gi = g.id AND src.strain = str.id AND str.gbuid = ? AND src.start = ? AND src.end = ?;
                ''', (strain_gbuid, start, end)).fetchone(), True)

        if not gi:
            gis[(name, gi_gbuid)] += 1
            n = name
            if gis[(name, gi_gbuid)] > 1:  # Disambiguate duplicate GI names
                if name != "cGI":
                    n = name + f" (cGI{gis[(name, gbuid)]})"
                else:
                    n = name = f"cGI{gis[(name, gbuid)]}"
            vals = dict(name=n, type=type, role=role)
            try:
                gi = conn.execute(f'''INSERT INTO genomic_islands (name, type, role) VALUES (?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                # This is completely useless given the appended (cGI#) to the name
                gi = dup(line, 'genomic_islands', vals, conn.execute(f'''SELECT * FROM genomic_islands WHERE name = ?;''', (name,)).fetchone())

        if gi is None:
            print(line, 'No GI, skipping', row)
            continue

        conn.execute(f'''INSERT OR IGNORE INTO alternate_names (gi, name) VALUES (?,?);''', (gi, name))

        # strains
        strain = None
        name = row[4].strip().replace(' ', ' ')  # TODO grab strain name from ncbi based on accession
        strain_gbuid = row[5].strip() or None
        vals = dict(name=name, gbuid=strain_gbuid)
        if name:
            try:
                strain = conn.execute(f'''INSERT INTO strains (name, gbuid) VALUES (?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                strain = dup(line, 'strains', vals, conn.execute(f'''SELECT id, name, gbuid FROM strains WHERE name = ? AND gbuid = ? UNION SELECT id, name, gbuid FROM strains WHERE gbuid = ? LIMIT 1;''', (name, strain_gbuid, strain_gbuid)).fetchone())
        else:
            print(line, "No strain name, skipping strain", row)

        # sequence
        seq = None
        seq_gc = None
        try: gc = float(row[3].strip()) or None
        except: gc = None
        paths = row[13].replace(',', ';').replace('and', ';').split(';') or None
        paths = list(p for p in paths if p.strip())
        gbuid = row[6].strip() or None
        if gc and not paths and not gbuid:
            print(f"{line} GC provided but no path, dropping GC data")
        if not paths and gbuid:
            paths.append(gbuid + '.fna')
        for path in paths:
            path = path.strip() or None
            if path:
                path = 'sequences/' + path
                if not os.path.isfile(path) and gbuid:
                    seq_gc = fetch(line, gbuid, path)
                    if seq_gc is None:
                        path = None
                    elif gc is None:
                        gc = seq_gc
                elif os.path.isfile(path):
                    for record in SeqIO.parse(path, 'fasta'):
                        seq_gc = round(GC(record.seq), 1)  # should only be one record
            if gc is not None and gc != seq_gc:
                print(f"{line} Calculated gc {seq_gc} doesn't match stored gc {gc}")
            vals = dict(gi=gi, gbuid=gbuid, gc=gc, path=path)
            try:
                seq = conn.execute(f'''INSERT INTO sequences (gi, gbuid, gc, path) VALUES (?,?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                seq = dup(line, 'sequences', vals, conn.execute(f'''SELECT id, gi, gbuid, gc, path FROM sequences WHERE gbuid = ?;''', (gbuid,)).fetchone())

        # source
        source = None
        try: size = int(float(row[7]) * 1000) if row[7] else None
        except: size = None
        try: start = int(row[8]) if row[8] else None
        except: start = None
        try: end = int(row[9]) if row[9] else None
        except: end = None
        if not seq and strain_gbuid and start is not None and end is not None:
            id = '_'.join((strain_gbuid, str(start), str(end)))
            path = 'sequences/' + id + '.fna'
            if not os.path.isfile(path):
                gc = fetch(line, strain_gbuid, path, slice(start, end + 1), id=id, description=role)
                if gc is None:
                    print(line, 'Unable to fetch strain for source', strain_gbuid)
                    path = None
            if path:
                vals = dict(gi=gi, gc=gc, path=path)
                try:
                    seq = conn.execute(f'''INSERT INTO sequences (gi, gc, path) VALUES (?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
                except sqlite3.IntegrityError:
                    seq = dup(line, 'source sequence', vals, conn.execute(f'''SELECT id, gi, gc, path FROM sequences WHERE path = ?;''', (path,)).fetchone())
        source = conn.execute(f'''INSERT INTO sources (gi, strain, size, start, end, seq) VALUES (?,?,?,?,?,?) RETURNING id;''', (gi, strain, size, start, end, seq)).fetchone()['id']

        # publications
        publications = row[10].replace(',', ';').replace('and', ';').split(';')
        for publication in publications:
            try:
                conn.execute(f'''INSERT INTO publications (source, publication) VALUES (?, ?);''', (source, publication.replace(' ', ' ').strip()))
            except sqlite3.IntegrityError:
                pass  # not much to do in this case

        # pmids
        pmids = row[11].replace(',', ';').replace('and', ';').split(';')
        if re.search('^\d+$', row[12]):
            pmids.append(row[12])
        for pmid in pmids:
            pmid = pmid.strip()
            if pmid:
                try:
                    conn.execute(f'''INSERT INTO pmids (source, pmid) VALUES (?,?);''', (source, pmid))
                except sqlite3.IntegrityError:
                    pass  # not much to do in this case

        print(f'{line} gi:{gi} strain:{strain} seq:{seq} source:{source} associated')

    for (name, gbuid), count in gis.items():
        if count > 1:
            conn.execute(f'''UPDATE genomic_islands SET name = ? WHERE name = ?;''', (name + " (cGI1)", name))

conn.commit()
with open('strain_names.dump', mode='wb') as f:
    pickle.dump(strain_names, f)

print("Validating database...")
seqs = dict()
for id, gc, path, name, gi in conn.execute(f'''SELECT s.id, s.gc, s.path, g.name, s.gi from sequences s LEFT JOIN genomic_islands g on s.gi = g.id;'''):
    # TODO attempt to recover gbuid by searching sequence on NCBI
    if not os.path.isfile(path):
        print(f"sequence {id}: {path} does not exist")
        continue
    records = list(SeqIO.parse(path, 'fasta'))
    if len(records) != 1:
        print(f'')
    for record in records:
        # verify seq.id uniqueness and accuracy
        if record.id in seqs:
            print(f"sequence {id}: duplicate sequence accession ({record.id}) found in {path} and {seqs[record.id]}")
        else:
            seqs[record.id] = path
            # change accessions in fastas to sequence ids to backlink to database from blast output
            record.id = f'{name}|{gi}|{id}'

        seq_gc = round(GC(record.seq), 1)
        if gc != seq_gc:
            print(f"sequence {id}: Calculated gc {seq_gc} doesn't match stored gc {gc}")

    with open(path, 'w') as f:
        SeqIO.write(records, f, 'fasta')

