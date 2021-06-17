import csv
import sqlite3
import re
import os

os.remove('GIs.sqlite')
conn = sqlite3.connect('GIs.sqlite')
conn.row_factory = sqlite3.Row
with open('schema.sql') as f:
    conn.executescript(f.read())

with open('GIs.csv') as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        # genomic_islands
        name = row[0].strip()
        if not name:
            print("No GI name")
            print(row)
            continue
        gi = conn.execute(f'''INSERT INTO genomic_islands (name, type, role) VALUES (?,?,?) RETURNING id;''', (name, row[1].strip() or None, row[2].strip() or None)).fetchone()['id']

        # strains
        name = row[4].strip()
        if not name:
            print("No strain name")
            print(row)
            continue
        strain = conn.execute(f'''INSERT INTO strains (name, gbuid) VALUES (?,?) RETURNING id;''', (name, row[5].strip() or None)).fetchone()['id']

        # sequence
        seq = None
        gc = row[3].strip() or None
        path = row[13].strip() or None
        if bool(gc) != bool(path):
            print("GC provided but no path")
            print(row)
            continue
        if gc and path:
            seq = conn.execute(f'''INSERT INTO sequences (gi, gbuid, gc, path) VALUES (?,?,?,?) RETURNING id;''', (gi, row[6].strip() or None, gc, path)).fetchone()['id']

        # source
        size = int(row[7]) if row[7] else None
        start = int(row[8]) if row[8] else None
        end = int(row[9]) if row[9] else None
        source = conn.execute(f'''INSERT INTO sources (gi, strain, size, start, end, seq, publication) VALUES (?,?,?,?,?,?,?) RETURNING id;''', (gi, strain, size, start, end, seq, row[10].strip() or None)).fetchone()['id']

        # pmids
        pmids = row[11].split(',')
        pmids.extend(pmids[0].split('and'))
        pmids[0] = ''
        if re.search('^\d+$', row[12]):
            pmids.append(row[12])
        for pmid in pmids:
            pmid = pmid.strip()
            if pmid:
                conn.execute(f'''INSERT INTO pmids (source, pmid) VALUES (?,?);''', (source, pmid))



