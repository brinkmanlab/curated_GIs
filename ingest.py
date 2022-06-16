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
Entrez.api_key = '3aa8b3eeefd28c87a8894c592cf54f305f09'

stdout = sys.stdout
sys.stdout = open('ingest.log', 'w')

strain_names = {}
if os.path.isfile('strain_names.dump'):
    with open('strain_names.dump', 'rb') as f:
        strain_names = pickle.load(f)

clean = False  # Delete database and start clean during import
fileext = 'gb'  # or 'fasta'


def fetch(line, gbuid, path, slice=None, id=None, description=None):
    """
    Fetch sequence associated with gbuid from NCBI nucleotide database
    :param line: CSV record line #
    :param gbuid: gbuid of sequence to fetch
    :param path: path to write fasta
    :param slice: slice to apply to sequence before writing
    :param id: id to assign to sequence in fasta
    :param description: description to include with sequence in fasta
    :return: gc content of sequence
    """
    try:
        handle = Entrez.efetch(db='nucleotide', rettype='gbwithparts' if fileext == 'gb' else fileext, retmode='text', id=gbuid)
        records = list(SeqIO.parse(handle, fileext))
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
            length = len(records[0])
            if slice.start >= length or slice.stop > length or slice.stop - slice.start > length:
                print(line, gbuid, "Warning: Sequence returned by NCBI is smaller than provided coordinates", length, slice)
            if slice.start >= slice.stop:
                # Crosses origin
                records[0] = records[0][slice.start:] + records[0][:slice.stop]
            else:
                records[0] = records[0][slice]
        if id is not None:
            records[0].id = id
        if description is not None:
            records[0].description = description
        # scrub bad chars from description
        records[0].description = records[0].description.replace('ﬄ', 'ffl').replace('\u2010', '-')
        with open(path, 'w') as file:
            SeqIO.write(records, file, fileext)
        return round(GC(records[0].seq), 1), len(records[0].seq)
    except HTTPError as e:
        print(line, gbuid, "No sequence returned by NCBI or bad request", e)
        return None


def get_strain_name(line, gbuid, name):
    """
    Fetch the strain name as listed in the NCBI taxonomy database
    :param line: CSV record line #
    :param gbuid: gbuid of strain
    :param name: name of strain in CSV
    :return: NCBI name, or original name if not found
    """
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
            taxid = int(taxid)
            records = Entrez.read(handle)
            handle.close()
            if len(records) < 1:
                print(line, gbuid, "Strain taxid returned no results by NCBI", taxid)
                return name
            if len(records) > 1:
                print(line, gbuid, "More than one record returned by NCBI for taxid", taxid, tuple(r.id for r in records))
                return name
            if records[0]['Rank'] != 'strain':
                print(line, gbuid, "Accession associated with non-strain taxa", taxid, f"'{records[0]['Rank']}', keeping", name)
                return name
            new_name = str(records[0]['ScientificName'])
            strain_names[gbuid] = new_name
        except (HTTPError, RuntimeError) as e:
            print(line, gbuid, "No records returned by NCBI or bad request:", e)
            return name
    if name != new_name:
        print(line, gbuid, f"Provided strain name ({name}) replaced with NCBI strain name ({records[0]['ScientificName']})")
    return new_name


def dup(line, table, new, old, quiet=False):
    """Handle duplicate unique columns and other table contraints
    :param line: CSV Line #
    :param table: Table name
    :param new: Row values being inserted
    :param old: Row values already in table
    :param quiet: Do not output non-duplcate errors
    :return: id generated from insertion
    """
    if old is None:
        if not quiet:
            print(f'{line} {table}: failed to meet constraint, discarding {tuple(new.values())}')
        return None
    elif any(old[i] != new[i] for i in new.keys()):
        print(f'{line} {table}: duplicate  {tuple(old or (old,))}')
        print(f'{line} {table}: discarding {tuple(new.values())}')
    return old['id']


def loadcsv(reader, conn):
    """
    Load CSV into database schema
    :param reader: CSV reader iterator
    :param conn: sqlite database connection
    :return: None
    """
    next(reader)  # skip header
    line = 1
    gis = defaultdict(int)
    for row in reader:
        names, type, role, _, strain_name, strain_gbuid, gi_gbuid, size, start, end, publications, _, _, _, *row_remainder = row
        line += 1

        # Handle coordinates with ':'
        start = start.split(':')
        strain_start = start[0]
        gi_start = None
        if len(start) > 1:
            gi_start = start[1]
            try: gi_start = int(gi_start) if gi_start else None
            except: gi_start = None
        try: strain_start = int(strain_start) if strain_start else None
        except: strain_start = None
        end = end.split(':')
        strain_end = end[0]
        gi_end = None
        if len(end) > 1:
            gi_end = end[1]
            try: gi_end = int(gi_end) if gi_end else None
            except: gi_end = None
        try: strain_end = int(strain_end) if strain_end else None
        except: strain_end = None
        del start, end
        gi_gbuid = gi_gbuid.strip() or None
        strain_gbuid = strain_gbuid.strip() or None
        strain_seq_id = None
        path = None
        if gi_gbuid:
            path = f"sequences/{gi_gbuid}.{fileext}"
        elif strain_gbuid is not None and strain_start is not None and strain_end is not None:
            strain_seq_id = '_'.join((strain_gbuid, str(strain_start), str(strain_end)))
            path = 'sequences/' + strain_seq_id + f".{fileext}"
        try: size = int(float(size) * 1000) if size else None
        except: size = None

        # genomic_islands
        gi = None
        names = names.strip().replace(' ', ' ').split('/')
        name = names[0]
        if not name:
            print(line, "No GI name, substituting 'cGI'", row)
            name = "cGI"
            names = [name]
        type = type.strip().replace(' ', ' ').lower().replace('-','_').replace(' ', '_') or None
        if type in ('bacteriophage', 'phage', 'prophage_like', 'putative_prophage'): type = 'prophage'
        if type not in ("prophage","ICE","transposon","putative_prophage","prophage_like","integron","integrated_plasmid",None):
            print(f"{line} discarding gi type {type}")
            type = None
        role = role.strip().replace(' ', ' ') or None
        vals = dict(name=name, type=type, role=role)
        # deduplicate gi on source accession and start-end or gi accession
        if gi_gbuid:
            gi = dup(line, 'genomic_islands', vals, conn.execute(f'''SELECT genomic_islands.* FROM genomic_islands, sequences WHERE gbuid = ? AND genomic_islands.id = sequences.gi''', (gi_gbuid,)).fetchone(), True)
        if not gi and strain_gbuid and strain_start is not None and strain_end is not None:
            gi = dup(line, 'genomic_islands', vals, conn.execute(f'''
                SELECT g.* FROM genomic_islands as g, sources as src, sequences as seq, strains as str 
                WHERE src.gi = g.id AND src.strain = str.id AND str.gbuid = ? AND src.start = ? AND src.end = ?;
                ''', (strain_gbuid, strain_start, strain_end)).fetchone(), True)

        if not gi:
            gis[name] += 1
            n = name
            if gis[name] > 1:  # Disambiguate duplicate GI names
                if name != "cGI":
                    n = name + f" (cGI{gis[name]})"
                else:
                    n = name = f"cGI{gis[name]}"
            vals = dict(name=n, role=role)
            try:
                gi = conn.execute(f'''INSERT INTO genomic_islands (name, role) VALUES (?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                gi = None #dup(line, 'genomic_islands', vals, conn.execute(f'''SELECT * FROM genomic_islands WHERE name = ?;''', (name,)).fetchone())
                print(line, 'GI name already exists, failed to increment cGI suffix')

        if gi is None:
            print(line, 'No GI, skipping', row)
            continue

        for n in names:
            conn.execute(f'''INSERT OR IGNORE INTO alternate_names (gi, name) VALUES (?,?);''', (gi, n))

        if type:
            conn.execute(f'''INSERT OR IGNORE INTO gi_type (gi, type) VALUES (?,?);''', (gi, type))

        # strains
        strain = None
        strain_name = strain_name.strip().replace(' ', ' ')
        strain_name = get_strain_name(line, strain_gbuid, strain_name)
        vals = dict(name=strain_name, gbuid=strain_gbuid)
        if strain_name:
            try:
                strain = conn.execute(f'''INSERT INTO strains (name, gbuid) VALUES (?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError as e:
                strain = dup(line, 'strains', vals, conn.execute(f'''SELECT id, name, gbuid FROM strains WHERE name = ? AND gbuid = ? UNION SELECT id, name, gbuid FROM strains WHERE gbuid = ? UNION SELECT id, name, gbuid FROM strains WHERE name = ? AND gbuid IS NULL LIMIT 1;''', (strain_name, strain_gbuid, strain_gbuid, strain_name)).fetchone())
        else:
            print(line, "No strain name, skipping strain", row)

        # sequence
        seq = None
        seq_gc = None
        seq_len = size
        seq_slice = None
        if path:
            if not os.path.isfile(path):
                print(f"{line} Missing {path}, attempting to download..")
                if gi_gbuid:
                    if gi_start is not None and gi_end is not None:
                        seq_slice = slice(gi_start-1, gi_end)
                    seq_gc, seq_len = fetch(line, gi_gbuid, path, seq_slice)
                elif strain_seq_id:
                    if strain_start is not None and strain_end is not None:
                        seq_slice = slice(strain_start - 1, strain_end)
                    seq_gc, seq_len = fetch(line, strain_gbuid, path, seq_slice, id=strain_seq_id, description=role)
                if seq_gc is None:
                    print(line, "Unable to download", path)
                    path = None
            else:
                for record in SeqIO.parse(path, fileext):
                    seq_gc = round(GC(record.seq), 1)  # should only be one record
                    seq_len = len(record.seq)
        if gi_gbuid or seq_gc or seq_len or path:
            vals = dict(gi=gi, gbuid=gi_gbuid, start=gi_start, end=gi_end, gc=seq_gc, length=seq_len, path=path)
            try:
                seq = conn.execute(f'''INSERT INTO sequences (gi, gbuid, start, end, gc, length, path) VALUES (?,?,?,?,?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                seq = dup(line, 'sequences', vals, conn.execute(f'''SELECT id, gi, gbuid, start, end, gc, length, path FROM sequences WHERE gbuid = ?;''', (gi_gbuid,)).fetchone())


        # source
        source = None
        vals = dict(gi=gi, strain=strain, start=strain_start, end=strain_end, seq=seq)
        try:
            source = conn.execute(f'''INSERT INTO sources (gi, strain, start, end, seq) VALUES (?,?,?,?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
        except sqlite3.IntegrityError:
            source = dup(line, 'sources', vals, conn.execute(f'''SELECT id, gi, strain, start, end, seq FROM sources WHERE strain = ? and start = ? and end = ?;''', (strain, strain_start, strain_end)).fetchone())

        # publications
        publications = publications.split(';')  # .replace(',', ';')
        for publication in publications:
            publication = publication.replace(' ', ' ').strip()
            if not publication:
                continue
            doi = (re.search('(?<=doi\.org/)\S+|(?<=doi:)\S+|(?<=doi: )\S+', publication) or [''])[0].rstrip('.') or None
            vals = dict(publication=publication, doi=doi)
            try:
                publication_id = conn.execute(f'''INSERT INTO publications (publication, doi) VALUES (?,?) RETURNING id;''', tuple(vals.values())).fetchone()['id']
            except sqlite3.IntegrityError:
                publication_id = dup(line, 'publications', vals, conn.execute(f'''SELECT id, publication, doi FROM publications WHERE publication = ? OR doi = ?;''', tuple(vals.values())).fetchone())
            conn.execute(f'''INSERT INTO source_pub_assoc (source, publication) VALUES (?, ?);''', (source, publication_id))

        print(f'{line} gi:{gi} strain:{strain} seq:{seq} source:{source} associated')

    for name, count in gis.items():
        if count > 1:
            conn.execute(f'''UPDATE genomic_islands SET name = ? WHERE name = ?;''', (name + " (cGI1)", name))

    conn.commit()
    with open('strain_names.dump', mode='wb') as f:
        pickle.dump(strain_names, f)


def validate(conn):
    """
    Scan through database, validating various aspects of the data
    Also changes sequence id to "gi_name|gi_id|seq_id" in file
    :param conn: sqlite database connection
    :return: None
    """
    print("Validating database...")
    seqs = dict()
    paths = []
    for id, gc, path, name, gi, strain, gbuid, strain_gbuid, start, end, seq_start, seq_end in conn.execute(f'''SELECT s.id, s.gc, s.path, g.name, s.gi, s3.name, s.gbuid, s3.gbuid, s2.start, s2.end, s.start, s.end from sequences s LEFT JOIN genomic_islands g on s.gi = g.id left join sources s2 on s.id = s2.seq left join strains s3 on s3.id = s2.strain;'''):
        if path is None:
            print(f"sequence {id}: has no associated path")
            continue
        paths.append(path)
        if not os.path.isfile(path):
            print(f"sequence {id}: {path} does not exist")
            continue
        records = list(SeqIO.parse(path, fileext))
        if len(records) != 1:
            print(f'multiple records found in {path}')
        for record in records:
            # change accessions in fastas to sequence ids to backlink to database from blast output
            record.id = record.name = f'{name.replace(" ", "_")}|{gi}|{id}'
            record.description = strain or 'No associated strain'
            # verify seq.name uniqueness and accuracy
            if record.name in seqs:
                print(f"sequence {id}: duplicate sequence accession ({record.name}) found in {path} and ({seqs[record.name][0]}) {seqs[record.name][1]}")
            seqs[record.name] = [id, path]

            seq_gc = round(GC(record.seq), 1)
            if gc != seq_gc:
                print(f"sequence {id}: Calculated gc {seq_gc} doesn't match stored gc {gc}")

            if len(record.seq) > 1000000:
                print(f"Warning: {path} contains a sequence greater than 1Mb in length")

            # check provided file if they match a downloaded copy
            #if gbuid or (strain_gbuid and start is not None and end is not None):
            #    if not os.path.isfile(path+'.tmp'):
            #        if gbuid:
            #            fetch(id, gbuid, path + '.tmp', slice(seq_start-1, seq_end) if seq_start is not None and seq_end is not None else None)
            #        else:
            #            fetch(id, strain_gbuid, path + '.tmp', slice(start-1, end))
            #    if os.path.isfile(path + '.tmp'):
            #        tmp_record = next(SeqIO.parse(path+'.tmp', fileext))
            #        if tmp_record.seq != record.seq:
            #            print(f"Warning: Sequence mismatch between NCBI and {path} ({len(tmp_record.seq)}bp vs {len(record.seq)}bp)")
            #elif not (gbuid or (strain_gbuid and start is not None and end is not None)):
            #    print(f"No strain or GI accession, unable to download comparison for {path}")

            #if not gbuid:
            #    # TODO attempt to recover gbuid by searching sequence on NCBI
            #    pass

        with open(path, 'w') as f:
            SeqIO.write(records, f, fileext)

    extra = set(f'sequences/{f}' for f in os.listdir('sequences/') if f.endswith(fileext)) - set(paths)
    if len(extra) > 0:
        extra = '\n'.join(extra)
        print(f"files detected in sequences/ that are not present in the database:\n{extra}")


if __name__ == "__main__":
    if clean:
        os.remove('GIs.sqlite')
    conn = sqlite3.connect('GIs.sqlite')
    conn.row_factory = sqlite3.Row
    if clean:
        with open('schema.sql') as f:
            conn.executescript(f.read())
        with open('GIs.csv') as f:
            reader = csv.reader(f)
            loadcsv(reader, conn)
    validate(conn)