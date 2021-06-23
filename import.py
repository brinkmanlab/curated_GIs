import csv
import sqlite3
import re
import os
import sys

stdout = sys.stdout
sys.stdout = open('ingest.log', 'w')
os.remove('GIs.sqlite')
conn = sqlite3.connect('GIs.sqlite')
conn.row_factory = sqlite3.Row
with open('schema.sql') as f:
    conn.executescript(f.read())


def dup(new, old):
    if any(old[i] != new[i] for i in new.keys()):
        print(f'duplicate  {tuple(old)}')
        print(f'discarding {tuple(new.values())}')
    return old['id']


with open('GIs.csv') as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        # genomic_islands
        name = row[0].strip().replace(' ', ' ')
        type = row[1].strip().replace(' ', ' ') or None
        role = row[2].strip().replace(' ', ' ') or None
        vals = dict(name=name, type=type, role=role)
        if not name:
            print("No GI name")
            print(row)
            continue
        try:
            gi = conn.execute(f'''INSERT INTO genomic_islands (name, type, role) VALUES (?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
        except sqlite3.IntegrityError:
            gi = dup(vals, conn.execute(f'''SELECT * FROM genomic_islands WHERE name = ?;''', (name,)).fetchone())

        # strains
        name = row[4].strip().replace(' ', ' ')
        gbuid = row[5].strip() or None
        vals = dict(name=name, gbuid=gbuid)
        if name:
            try:
                strain = conn.execute(f'''INSERT INTO strains (name, gbuid) VALUES (?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                strain = dup(vals, conn.execute(f'''SELECT id, name, gbuid FROM strains WHERE name = ? OR gbuid = ?;''', tuple(vals.values())).fetchone())
        else:
            print("No strain")
            print(row)

        # sequence
        seq = None
        try: gc = float(row[3].strip()) or None
        except: gc = None
        paths = row[13].replace(',', ';').replace('and', ';').split(';') or None
        if len(paths) > 1:
            gc = None
        if bool(gc) != bool(paths) and False:
            print("GC provided but no path, dropping GC data")
            print(row)
        gbuid = row[6].strip() or None
        for path in paths:
            path = path.strip()
            if path:
                vals = dict(gi=gi, gbuid=gbuid, gc=gc, path='sequences/' + path)
                try:
                    seq = conn.execute(f'''INSERT INTO sequences (gi, gbuid, gc, path) VALUES (?,?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
                except sqlite3.IntegrityError:
                    seq = dup(vals, conn.execute(f'''SELECT id FROM sequences WHERE gbuid = ?;''', (gbuid,)).fetchone())

        # source
        try: size = int(float(row[7]) * 1000) if row[7] else None
        except: size = None
        try: start = int(row[8]) if row[8] else None
        except: start = None
        try: end = int(row[9]) if row[9] else None
        except: end = None
        source = conn.execute(f'''INSERT INTO sources (gi, strain, size, start, end, seq, publication) VALUES (?,?,?,?,?,?,?) RETURNING id;''', (gi, strain, size, start, end, seq, row[10].strip() or None)).fetchone()['id']

        # pmids
        pmids = row[11].replace(',', ';').replace('and', ';').split(';')
        if re.search('^\d+$', row[12]):
            pmids.append(row[12])
        for pmid in pmids:
            pmid = pmid.strip()
            if pmid:
                conn.execute(f'''INSERT INTO pmids (source, pmid) VALUES (?,?);''', (source, pmid))

conn.commit()

for id, path in conn.execute(f'''SELECT id, path from sequences;'''):
    if not os.path.isfile(path):
        print(f"sequence {id} {path} does not exist")