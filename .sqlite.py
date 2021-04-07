#!/usr/bin/env python3
"""
SQLite External diff and merge tool for Git

Execute this script within a repository to register it with Git
"""
import subprocess
import sys
import os
import sqlite3


def deltas(n, o, pk, existing_rows=True):
    """
    Build query for rows that differ between 'main' and 'old' schema for table
    :param n: new table fqn
    :param o: old table fqn
    :param pk: List of primary key names or None for natural join
    :param existing_rows: Set to false to include all rows during natural join
    :return: SQL query as string with
    """
    if pk is None:
        return f'''
        SELECT n.*, o.*
        FROM {n} as n
        NATURAL JOIN {o} as o
        ''' + ';' if existing_rows else f'''
        EXCEPT (
            SELECT n.*, o.*
            FROM {n} as n
            LEFT JOIN {o} as o WHERE 0
            UNION ALL
            SELECT n.*, o.*
            FROM {o} as o
            LEFT JOIN {n} as n WHERE 0
        );
        '''
    return f'''
        SELECT n.*, o.*
        FROM {n} as n
        INNER JOIN (
          SELECT {','.join(f'{n}.{k}' for k in pk)} FROM {n} NATURAL LEFT JOIN {o} WHERE {keys_null(o, pk)}
        ) as pk ON ({' AND '.join(f'n.{k} = pk.{k}' for k in pk)})
        INNER JOIN {o} as o ON ({' AND '.join(f'n.{k} = o.{k}' for k in pk)});
    '''


def attach(conn, databases):
    """
    Attach a mapping of scheme names to file paths to the connection
    :param conn: SQLite connection
    :param databases: dict of paths to sqlite databases, keyed on scheme name
    :return: None
    """
    for scheme, db in databases.items():
        conn.execute(f"ATTACH '{db}' AS {scheme};")


def get_tables(conn, scheme):
    """
    Gather table data into a dict of Table objects for a specific scheme in the connection
    :param conn: SQLite database connection
    :param scheme: name of scheme to read
    :return: dict of Table objects keyed on name
    """
    result = {}
    idx = []
    cur = conn.execute(f'SELECT * FROM {scheme}.sqlite_master;')
    for table in cur:
        if table['name'].startswith('sqlite_'):
            continue
        if table['type'] == 'table':
            result[table['name']] = Table(conn, scheme, table['name'], table['sql'])
        elif table['type'] == 'index' and table['sql']:
            idx.append((table['name'], table['tbl_name'], table['sql']))
        elif table['type'] == 'view':
            result[table['name']] = View(conn, scheme, table['name'], table['sql'])
    for name, table, sql in idx:
        result[table].indexes[name] = sql
    return result


def to_sql(val):
    """
    Convert a python native value to SQL
    :param val: value to convert
    :return: string of SQL value ready to insert in query
    """
    if val is None:
        return 'NULL'
    return ("'" + val.replace("'", "''") + "'") if isinstance(val, str) else val


def keys_null(table, keys):
    """
    Generate a SQL condition statement checking if all keys are null
    :param table: table the statement applies to
    :param keys: list of column names to check if null
    :return: string containing SQL conditional
    """
    return ' AND '.join(f'{table}.{k} IS NULL' for k in keys)


def USING(n, o, pk):
    """
    SQLite does not work with USING(rowid), this function replaces that
    :param n: name of table1
    :param o: name of table2
    :param pk: list of primary keys to use in join
    :return: string containing ON clause equivalent to USING
    """
    return f'''ON ({' AND '.join(f'{n}.{k}={o}.{k}' for k in pk)})'''


def prompt(message, values):
    """
    Prompt user for single letter choice
    :param message: Message to display to user
    :param values: valid choices
    :return: users choice
    """
    while True:
        val = input(message)
        if val[0] in values:
            return val[0]
        else:
            print(f"Valid options: {', '.join(values)}")


def dump(dump_path, working_sql, ancestor, working, remote, modified_tables):
    """
    Write SQL file with header to be read by resume()
    :param dump_path: Path to write output to
    :param working_sql: path to existing sql file
    :param ancestor: path to ancestor sqlite database
    :param working: path to working sqlite database
    :param remote: path to remote sqlite database
    :param modified_tables: set of table names that are included in the working_sql
    :return: None
    """
    print('Edit, resolve conflicts, and run', dump_path)
    with open(dump_path, 'w') as dump_file, open(working_sql, 'r') as working_schema:
        print(f'#!{sys.argv[0]} resume {ancestor} {working} {remote}', file=dump_file)
        print('"' + '","'.join(modified_tables) + '"')
        dump_file.writelines(working_schema)
        os.chmod(dump_file.fileno(), 0o774)


class Column:
    def __init__(self, cid, name, type, notnull, default, pk):
        self.cid = cid
        self.name = name
        self.type = type
        self.notnull = bool(notnull)
        self.default = default
        self.pk = bool(pk)

    def __eq__(self, other):
        if isinstance(other, str):
            return self.name == other
        if not isinstance(other, Column):
            return False
        return self.cid == other.cid and self.name == other.name and self.type == other.type and self.notnull == other.notnull and self.default == other.default and self.pk == other.pk

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return f'"{self.name}"'


class Table:
    def __init__(self, conn, scheme, name, sql):
        self.scheme = scheme
        self.name = name
        self.fqn = f'"{self.scheme}"."{self.name}"'
        self.sql = sql
        cur = conn.execute(f"PRAGMA \"{scheme}\".TABLE_INFO('{name}')")
        self.pk = []
        self.cols = {}
        self.indexes = {}
        for col in cur:
            c = Column(*col)
            self.cols[c.name] = c
            if c.pk:
                self.pk.append(c)
        if not len(self.pk):
            self.pk = ('rowid',)
        else:
            self.pk = tuple(self.pk)
        self.hasRowID = 'without rowid' not in self.sql.lower()
        self.compositePk = len(self.pk) > 1

    def __str__(self):
        return f'"{self.scheme}"."{self.name}"'

    def __eq__(self, other):
        if not isinstance(other, Table):
            return False
        if len(self.cols) != len(other.cols):
            return False
        # Compare SQL to catch comment changes
        # for name, col in self.cols.items():
        #     if name not in other.cols:
        #         return False
        #     if col != other.cols[name]:
        #         return False
        return self.sql == other.sql

    def __ne__(self, other):
        return not self.__eq__(other)

    def dump_rows(self, rows):
        col_names = ','.join(f'"{name}"' for name in self.cols.keys())
        for row in rows:
            vals = ','.join(to_sql(row[name]) for name in self.cols.keys())
            print(f'INSERT INTO "{self.name}"({col_names}) VALUES ({vals});')

    def compare_indexs(self, other):
        new_indexes = set(self.indexes.keys())
        old_indexes = set(other.keys())
        created_indexes = new_indexes - old_indexes
        dropped_indexes = old_indexes - new_indexes
        for name in new_indexes.intersection(old_indexes):  # Drop and create modified indexes
            if self.indexes[name] != other[name]:
                created_indexes.add(name)
                dropped_indexes.add(name)
        return created_indexes, dropped_indexes


class View:
    def __init__(self, conn, scheme, name, sql):
        self.scheme = scheme
        self.name = name
        self.sql = sql
        self.indexes = {}

    def __str__(self):
        return f'"{self.scheme}"."{self.name}"'

    def __eq__(self, other):
        if not isinstance(other, View):
            return False
        return self.sql == other.sql

    def __ne__(self, other):
        return not self.__eq__(other)


def diff(path, old_file, old_hex, old_mode, new_file, new_hex, new_mode):
    conn = sqlite3.connect(new_file)
    conn.row_factory = sqlite3.Row
    attach(conn, {"old": old_file})
    new_tables = get_tables(conn, 'main')
    old_tables = get_tables(conn, 'old')

    new_table_names = set(new_tables.keys())
    old_table_names = set(old_tables.keys())

    created_tables = new_table_names - old_table_names
    dropped_tables = old_table_names - new_table_names

    # Print header
    print(f".sqlite.py diff a/{os.path.basename(old_file)} b/{os.path.basename(new_file)}")
    mode = ""
    if old_mode != new_mode:
        print(f"old mode {old_mode}")
        print(f"new mode {new_mode}")
    else:
        mode = f" {new_mode}"

    print(f"index {old_hex[0:7]}..{new_hex[0:7]}{mode}")

    # Print dropped tables
    for table in dropped_tables:
        print(f'DROP TABLE "{table}";')

    # Check if existing tables have been modified
    existing_tables = new_table_names.intersection(old_table_names)
    for table in existing_tables:
        new_table_type = 'TABLE' if isinstance(new_tables[table], Table) else 'VIEW'
        old_table_type = 'TABLE' if isinstance(old_tables[table], Table) else 'VIEW'
        if new_tables[table] != old_tables[table]:  # table/view modified
            print(f'DROP {old_table_type} "{table}"; -- schema change')
            print(new_tables[table].sql)

            if new_table_type == 'TABLE':
                #dump all rows, noting row deltas
                pk = new_tables[table].pk
                if new_tables[table].pk != old_tables[table].pk:
                    if new_tables[table].hasRowID and old_tables[table].hasRowID:
                        pk = ('rowid',)
                    else:
                        # Natural join, because can't compare PKs
                        pk = None

                cur = conn.execute(deltas(new_tables[table].fqn, old_tables[table].fqn, pk, False))
                for row in cur:
                    keys = row.keys()[:len(new_tables[table].cols)]
                    col_names = ",".join(f'"{key}"' if key.lower() != 'rowid' else key for key in keys)
                    new_vals = row[:len(new_tables[table].cols)]
                    old_vals = row[len(new_tables[table].cols):]
                    vals = ",".join(to_sql(val) for val in new_vals)
                    if all(v is None for v in new_vals):
                        # Row dropped
                        print(f'''-- DELETE FROM "{table}" WHERE {",".join(f'{k} = {to_sql(v)}' for k, v in zip(keys, new_vals))}''')
                    elif all(v is None for v in old_vals):
                        # Row added
                        print(f'INSERT INTO "{table}"({col_names}) VALUES ({vals});')
                    else:
                        # Row modified
                        old_vals = ",".join(f'{key}: {to_sql(val)}' for key, val in zip(row.keys()[len(new_tables[table].cols):], old_vals))
                        print(f'INSERT INTO "{table}"({col_names}) VALUES ({vals}); -- {{{old_vals}}}')
        elif old_table_type == 'TABLE':  # Don't diff VIEWs
            pk = old_tables[table].pk
            # new rows
            cur = conn.execute(f'''SELECT n.* FROM main."{table}" AS n LEFT JOIN old."{table}" AS o ON ({' AND '.join(f'n.{k}=o.{k}' for k in pk)}) WHERE {keys_null('o', pk)};''')
            new_tables[table].dump_rows(cur)

            # deleted rows
            cur = conn.execute(f'''SELECT o.* FROM old."{table}" AS o LEFT JOIN main."{table}" AS n ON ({' AND '.join(f'n.{k}=o.{k}' for k in pk)}) WHERE {keys_null('n', pk)};''')
            for row in cur:
                print(f'''DELETE FROM "{table}" WHERE {' AND '.join(f'{k} = {to_sql(row[k])}' for k in pk)};  -- ({", ".join(to_sql(v) for v in row)})''')

            # modified rows
            pk = new_tables[table].pk
            cur = conn.execute(deltas(new_tables[table].fqn, old_tables[table].fqn, pk))
            for row in cur:
                vals = []
                old_vals = []
                for i, (new_val, old_val) in enumerate(zip(row[:len(row)/2], row[len(row)/2:])):
                    if new_val != old_val:
                        vals.append(f'{row.keys()[i]} = {to_sql(new_val)}')  # assumes row.keys() is ordered same as table.*
                        old_vals.append(to_sql(old_val))
                print(f'''UPDATE "{table}" SET {", ".join(vals)} WHERE {' AND '.join(f'{k} = {to_sql(row[k])}' for k in pk)}; -- ({", ".join(old_vals)})''')

        # Compare indexes
        created_indexes, dropped_indexes = new_tables[table].compare_indexs(old_tables[table].indexes)
        for name in dropped_indexes:  # Drop deleted indexes
            print(f'DROP INDEX "{name}";')
        for name in created_indexes:  # Output created indexes
            print(new_tables[table].indexes[name])

    # Print new tables
    for name, table in new_tables.items():
        if name in created_tables:
            print(table.sql)
            for index in table.indexes.values():
                print(index)
            if not isinstance(table, View):
                cur = conn.execute(f'SELECT * FROM main."{name}"')
                new_tables[name].dump_rows(cur)
    return 0


def merge_table_data(conn, working_table, remote_table, ancestor_table, pk):
    if pk is None:
        # conflicts are unknowable, simply insert rows that are modified in remote
        cur = conn.execute(deltas(remote_table.fqn, ancestor_table.fqn, None))
        if cur.rowcount > 0 and prompt('Warning: Unable to match primary keys, some data may be duplicated in '+ working_table.name + '. Would you like to (c)ontinue or (s)kip?', ['c', 's']) == 'c':
            for row in cur:
                keys = row.keys()[:len(remote_table.cols)]
                col_names = ",".join(f'"{key}"' if key.lower() != 'rowid' else key for key in keys)
                new_vals = row[:len(remote_table.cols)]
                old_vals = row[len(remote_table.cols):]
                vals = ",".join(to_sql(val) for val in new_vals)
                if all(v is None for v in new_vals):
                    # Row dropped
                    pass  # Cant target row in working to drop
                elif all(v is None for v in old_vals):
                    # Row added
                    conn.execute(f'INSERT OR REPLACE INTO {working_table.fqn}({col_names}) VALUES ({vals});')
                else:
                    # Row modified
                    conn.execute(f'INSERT OR UPDATE INTO {working_table.fqn}({col_names}) VALUES ({vals});')
        return False
    PKS = (str(k) for k in pk)
    REMOTEPK = ",".join(f'r.{k}' for k in PKS)
    WORKPK = ",".join(f"w.{k}" for k in PKS)
    PK = ",".join(PKS)
    # Calc intersect of new or removed rows and emit conflicts
    # conflicts:
    # remote_new ^ working_new - only if different
    # remote_new ^ working_del
    # remote_new ^ working_mod
    # remote_del ^ working_new
    # remote_del ^ working_mod
    # remote_mod ^ working_new
    # remote_mod ^ working_del
    # remote_mod ^ working_mod - only if different
    conn.executescript(f'''
                DROP VIEW IF EXISTS temp.remote_new; 
                DROP VIEW IF EXISTS temp.working_new;
                DROP VIEW IF EXISTS temp.remote_del; 
                DROP VIEW IF EXISTS temp.working_del;
                DROP VIEW IF EXISTS temp.working_mod;
                DROP VIEW IF EXISTS temp.remote_mod; 
                DROP VIEW IF EXISTS temp.conflict;
                DROP TABLE IF EXISTS temp.working_conflicts;
                CREATE TEMP VIEW remote_new AS SELECT r.* FROM {remote_table.fqn} AS r LEFT JOIN {ancestor_table.fqn} AS a {USING('r','a', pk)} WHERE {keys_null('a', pk)};
                CREATE TEMP VIEW working_new AS SELECT w.* FROM {working_table.fqn} AS w LEFT JOIN {ancestor_table.fqn} AS a {USING('w','a', pk)} WHERE {keys_null('a', pk)};
                CREATE TEMP VIEW remote_del AS SELECT a.* FROM {ancestor_table.fqn} AS a LEFT JOIN {remote_table.fqn} AS r {USING('a','r', pk)} WHERE {keys_null('r', pk)};
                CREATE TEMP VIEW working_del AS SELECT a.* FROM {ancestor_table.fqn} AS a LEFT JOIN {working_table.fqn} AS w {USING('a','w', pk)} WHERE {keys_null('w', pk)};
                CREATE TEMP VIEW working_mod AS SELECT w.* FROM {working_table.fqn} AS w NATURAL LEFT JOIN {ancestor_table.fqn} AS a LEFT OUTER JOIN temp.working_new as p {USING('w','p', pk)} WHERE {keys_null('a', pk)} AND {keys_null('p', pk)};
                CREATE TEMP VIEW remote_mod AS SELECT r.* FROM {remote_table.fqn} AS r NATURAL LEFT JOIN {ancestor_table.fqn} AS a LEFT OUTER JOIN temp.remote_new as p {USING('r','p', pk)} WHERE {keys_null('a', pk)} AND {keys_null('p', pk)};

                CREATE TEMP VIEW conflict
                    AS    SELECT {REMOTEPK} FROM remote_new AS r INNER JOIN working_new AS w {USING('r','w',pk)}
                    EXCEPT SELECT {REMOTEPK} FROM remote_new AS r NATURAL INNER JOIN working_new AS w
                    UNION SELECT {REMOTEPK} FROM remote_new AS r INNER JOIN working_del AS w {USING('r','w',pk)}
                    UNION SELECT {REMOTEPK} FROM remote_new AS r INNER JOIN working_mod AS w {USING('r','w',pk)}
                    UNION SELECT {REMOTEPK} FROM remote_del AS r INNER JOIN working_new AS w {USING('r','w',pk)}
                    UNION SELECT {REMOTEPK} FROM remote_del AS r INNER JOIN working_mod AS w {USING('r','w',pk)}
                    UNION SELECT {REMOTEPK} FROM remote_mod AS r INNER JOIN working_new AS w {USING('r','w',pk)}
                    UNION SELECT {REMOTEPK} FROM remote_mod AS r INNER JOIN working_del AS w {USING('r','w',pk)}
                    UNION SELECT {REMOTEPK} FROM remote_mod AS r INNER JOIN working_mod AS w {USING('r','w',pk)}
                    EXCEPT SELECT {REMOTEPK} FROM remote_mod AS r NATURAL INNER JOIN working_mod AS w;
                ''')

    # new rows
    conn.execute(
        f'''INSERT OR IGNORE INTO {working_table.fqn} SELECT r.* FROM temp.remote_new AS r LEFT JOIN temp.conflict AS c {USING('r','c',pk)} WHERE {keys_null("c", pk)};''')

    # deleted rows
    conn.execute(f'DELETE FROM {working_table.fqn} WHERE ({PK}) IN (SELECT {PK} FROM temp.remote_del) AND ({PK}) NOT IN temp.conflict;')

    # modified rows without conflict
    cols = (f'"{col.name}"="{col.name}"' for col in working_table.cols)
    conn.execute(
        f'''UPDATE {working_table.fqn} AS w SET {", ".join(cols)} FROM temp.remote_mod WHERE ({WORKPK}) IS ({",".join(f"temp.remote_mod.{k}" for k in PKS)}) AND ({WORKPK}) NOT IN temp.conflict;''')

    # Resolve conflicting rows
    # create temporary copy of conflicts to allow modification while iterating
    conn.execute(f'''CREATE TEMP TABLE working_conflicts SELECT w.* FROM temp.conflict AS c INNER JOIN {working_table.fqn} AS w {USING('w','c',pk)};''')
    cur = conn.execute(f'''
                SELECT w.*, r.* FROM temp.working_conflicts AS w LEFT JOIN {remote_table.fqn} AS r {USING('w','r',pk)}
                UNION
                SELECT w.*, r.* FROM temp.conflict AS c INNER JOIN {remote_table.fqn} AS r {USING('r','c',pk)} LEFT JOIN temp.working_conflicts AS w {USING('c','w',pk)};
                ''')
    # Because the table columns are concatenated, index into query columns by original columns names
    working_table_col_i = {str(c): c.cid for c in working_table.cols.values()}
    remote_table_col_i = {str(c): c.cid for c in remote_table.cols.values()}
    split = len(working_table_col_i)
    for row in cur:  # iterate rows and resolve conflicts
        working_vals = {str(k): row[i] for k, i in working_table_col_i.items()}
        remote_vals = {str(k): row[i + split] for k, i in remote_table_col_i.items()}
        # deleted or created conflict
        working_missing = any(working_vals[str(k)] is None for k in pk)
        if working_missing or any(remote_vals[str(k)] is None for k in pk):  # Row missing from working or remote
            vals = working_vals if working_missing else remote_vals
            print((f'{k}: {to_sql(v)}' for k, v in vals.items()), sep=", ")
            action = prompt('(k)eep or (d)elete row? ', ('k', 'd'))
            if action == 'k' and working_missing:  # Only need to insert when missing from working
                keys = ",".join(vals.keys())
                values = ",".join(to_sql(v) for v in vals.values())
                conn.execute(f'INSERT OR IGNORE INTO {working_table.fqn}({keys}) VALUES ({values})')
            elif action == 'd' and working_missing:  # Only need to delete from working
                conn.execute(f'''DELETE FROM {working_table.fqn} ({PK}) = ({",".join("?" * len(pk))});''', (working_vals[str(k)] for k in pk))
        else:  # both rows exist, merge columns on name
            ancestor_vals = None
            unresolved_vals = {}
            final_vals = {}
            # Compare each row value to find differences and conflict
            for col in working_vals:
                if working_vals[col] != remote_vals[col]:
                    if ancestor_vals is None:  # Lazily fetch ancestor row
                        ancestor_vals = conn.execute(f'SELECT * FROM {ancestor_table.fqn} AS a WHERE ({PK}) = ({",".join("?" * len(pk))});',
                                                     (working_vals[str(k)] for k in pk)).fetchone() or False
                    if ancestor_vals:  # both working and remote updated
                        working_modified = working_vals[col] != ancestor_vals[col]
                        remote_modified = remote_vals[col] != ancestor_vals[col]
                        if working_modified and remote_modified:  # both modified
                            # Store unresolved value to later prompt user
                            unresolved_vals[col] = (working_vals[col], remote_vals[col], ancestor_vals[col], True)
                        elif remote_modified:  # only need to update if remote modified
                            final_vals[col] = remote_vals[col]
                    else:  # both new rows
                        # Store unresolved value to later prompt user
                        unresolved_vals[col] = (working_vals[col], remote_vals[col], None, False)
            # Prompt user to resolve conflict
            if unresolved_vals:
                print('Column value conflicts, choose the value to retain:')
                print('Working (', ', '.join(f'{k}: {to_sql(v)}' for k, v in working_vals.items()), ')')
                print('Remote (', ', '.join(f'{k}: {to_sql(v)}' for k, v in remote_vals.items()), ')')
                print('Ancestor (', ', '.join(f'{k}: {to_sql(v)}' for k, v in ancestor_vals.items()), ')')
                for col, (*v, has_ancestor) in unresolved_vals:
                    val_sep = ', '
                    last_sep = ''
                    if any(isinstance(o, str) and len(o) > 16 for o in v):
                        val_sep = last_sep = '\n'
                    w, r, a = v
                    ancestor_option = f'{val_sep}(a)ncestor: {a}' if has_ancestor else ''
                    action = prompt(f'(w)orking: {w}{val_sep}(r)emote: {r}{ancestor_option}{last_sep}{col}?',
                                    ('w', 'r', 'a') if has_ancestor else ('w', 'r'))
                    if action == 'w':
                        final_vals[col] = w
                    elif action == 'r':
                        final_vals[col] = r
                    elif action == 'a':
                        final_vals[col] = a
            # Update row with final values
            cols = (f'{k}={to_sql(v)}' for k, v in final_vals.items())
            conn.execute(f'''UPDATE {working_table.fqn} AS w SET {", ".join(cols)} WHERE {' AND '.join(f'w.{k} = {to_sql(working_vals[k])}' for k in pk)};''')
    return True


def merge(ancestor, working, remote, marker_size, placeholder):
    conn = sqlite3.connect(working)
    conn.row_factory = sqlite3.Row
    conn.isolation_level = 'EXCLUSIVE'
    attach(conn, {"ancestor": ancestor, "remote": remote})
    working_tables = get_tables(conn, 'main')
    remote_tables = get_tables(conn, 'remote')
    ancestor_tables = get_tables(conn, 'ancestor')

    working_table_names = set(working_tables.keys())
    remote_table_names = set(remote_tables.keys())

    # Tables in remote only
    remote_only_tables = remote_table_names - working_table_names
    for table in remote_only_tables:  # Copy table, indexes, and data into working
        conn.executescript(remote_tables[table].sql)
        for sql in remote_tables[table].indexes.values():  # copy indexes
            if sql:
                conn.executescript(sql)
        conn.execute(f'INSERT INTO main."{table}" FROM remote."{table}";')

    # Tables in both working and remote
    merge_tables = working_table_names.intersection(remote_table_names)
    modified_tables = set()  # Set of tables to merge sql with git

    # Build paths for working files to pass to git merge-file
    working_sql = os.path.basename(working) + '.working.sql'
    remote_sql = os.path.basename(remote) + '.remote.sql'
    ancestor_sql = os.path.basename(ancestor) + '.ancestor.sql'

    # Scan through tables and find tables that are modified
    modified_indexes = False
    for table in merge_tables:
        if working_tables[table] != remote_tables[table]:  # table modified
            # Resolve schema conflicts using sql to retain comments
            # Attempt to merge using git merge-file -L working -L ancestor -L remote --marker-size={marker_size} working ancestor remote
            for path, t in ((working_sql, working_tables[table]), (remote_sql, remote_tables[table]), (ancestor_sql, ancestor_tables[table])):
                with open(path, 'a') as schema:
                    schema.write(t.sql)
                    schema.write("\n")

            modified_tables.add(table)
        elif isinstance(working_tables[table], Table):  # No need to diff VIEWs
            merge_table_data(conn, working_tables[table], remote_tables[table], ancestor_tables[table], remote_tables[table].pk)
        # Compare indexes
        working_created_indexes, working_dropped_indexes = working_tables[table].compare_indexs(ancestor_tables[table].indexes)
        remote_created_indexes, remote_dropped_indexes = remote_tables[table].compare_indexs(ancestor_tables[table].indexes)
        conflicting_indexes = set(name for name in working_created_indexes.intersection(remote_created_indexes) if working_tables[table].indexes[name] != remote_tables[table].indexes[name])  # Find created or modified indexes that conflict
        for name in remote_dropped_indexes - conflicting_indexes:  # Drop deleted indexes
            conn.execute(f'DROP INDEX IF EXISTS "{name}";')
        for name in remote_created_indexes - conflicting_indexes:  # Output created indexes
            conn.execute(remote_tables[table].indexes[name])
        if conflicting_indexes:
            for path, t in ((working_sql, working_tables[table]), (remote_sql, remote_tables[table]), (ancestor_sql, ancestor_tables[table])):
                with open(path, 'a') as schema:
                    for name in conflicting_indexes:
                        schema.write(f'DROP INDEX IF EXISTS "{name}";\n')
                        schema.write(t.indexes[name])
                        schema.write("\n")
            modified_indexes = True

    # Merge modified table schemas with git merge-file
    if modified_tables or modified_indexes:
        dump_path = os.path.join(os.path.dirname(placeholder), os.path.basename(placeholder) + '.sql')
        merged = subprocess.run(["git", "merge-file", "-L", "working", "-L", "ancestor", "-L", "remote", "--marker-size", marker_size, working_sql, ancestor_sql, remote_sql], text=True, capture_output=True)
        if merged.returncode < 0:
            # git merge error
            print(merged.stderr, file=sys.stderr)
            dump(dump_path, working_sql, ancestor, working, remote, modified_tables)
        elif merged.returncode > 0:
            # merge conflicts
            git_editor = subprocess.run(["git", "config", "core.editor"], text=True, capture_output=True)
            git_editor_global = subprocess.run(["git", "config", "--global", "core.editor"], text=True, capture_output=True)
            editor = git_editor.stdout or git_editor_global.stdout or os.environ.get('VISUAL') or os.environ.get('EDITOR')
            if editor:
                action = prompt(f"Conflicts detected in table schemas. (e)dit or (w)rite to {dump_path}", ('e', 'w'))
                if action == 'e':
                    touch_time = os.stat(working_sql).st_mtime
                    subprocess.run([editor, working_sql], stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
                    if touch_time != os.stat(working_sql).st_mtime:
                        return resume(ancestor, working, remote, working_sql, modified_tables, conn, placeholder)
                    else:
                        print('No change detected')
                elif action == 'w':
                    dump(dump_path, working_sql, ancestor, working, remote, modified_tables)
            else:
                dump(dump_path, working_sql, ancestor, working, remote, modified_tables)
        else:
            # git merge successful
            return resume(ancestor, working, remote, working_sql, modified_tables, conn, placeholder)

        return merged.returncode
    return 0


def resume(ancestor, working, remote, dump_path, modified_tables=None, conn=None, placeholder=None):
    with open(dump_path, 'r') as dump_file:
        if modified_tables is None:
            dump_file.readline()  # skip shebang
            modified_tables = dump_file.readline()
            modified_tables = set(modified_tables[1:-1].split('","'))
        if conn is None:
            conn = sqlite3.connect(working)
            conn.row_factory = sqlite3.Row
            conn.isolation_level = 'EXCLUSIVE'
            attach(conn, {"ancestor": ancestor, "remote": remote})
        temp_prefix = "preparing_for_transform_"
        # Keep old tables to copy data
        incomplete = False
        try:
            for table in modified_tables:
                conn.execute(f'ALTER TABLE main."{table}" RENAME TO main."{temp_prefix}{table}";')
            conn.executescript(dump_file.read())
            working_tables = get_tables(conn, 'main')
            remote_tables = get_tables(conn, 'remote')
            ancestor_tables = get_tables(conn, 'ancestor')
            for table in modified_tables:
                # Copy in old working data into new schema
                new_cols = set(working_tables[table].cols.keys())
                old_cols = set(working_tables[temp_prefix+table].cols.keys())
                if new_cols == old_cols:
                    # Insert by column name
                    col_names = ",".join(f'"{c}"' for c in new_cols)
                    conn.execute(f'''INSERT INTO main."{table}" ({col_names}) SELECT {col_names} FROM main."{temp_prefix}{table}";''')
                elif len(new_cols) == len(old_cols):
                    # Insert by index
                    conn.execute(f'''INSERT INTO main."{table}" SELECT * FROM main."{temp_prefix}{table}";''')
                else:
                    # Insert subset of column names
                    common_columns = new_cols.intersection(old_cols)
                    col_names = ",".join(f'"{c}"' for c in common_columns)
                    conn.execute(f'''INSERT INTO main."{table}" ({col_names}) SELECT {col_names} FROM main."{temp_prefix}{table}";''')

                # Merge in remote data to new schema
                common_keys = set(working_tables[table].pk).intersection(set(working_tables[temp_prefix+table].pk))
                incomplete = not merge_table_data(conn, working_tables[table], remote_tables[table], ancestor_tables[table], common_keys or None) and incomplete
        except sqlite3.Error as e:
            print(e)
            if placeholder is not None:
                dump(os.path.join(os.path.dirname(placeholder), os.path.basename(placeholder) + '.sql'), dump_path, ancestor, working, remote, modified_tables)
            else:
                print(f"Please revise SQLite3 syntax and rerun {dump_path}")
            return 1
    return 1 if incomplete else 0


if __name__ == '__main__':
    if sys.argv[1] == 'diff':
        exit(diff(*(sys.argv[2:])))
    elif sys.argv[1] == 'merge':
        exit(merge(*(sys.argv[1:])))
    elif sys.argv[1] == 'resume':
        exit(resume(*(sys.argv[1:])))
    else:
        # Configure git
        subprocess.run(["git", "config", "diff.sqlite.binary", "true"])
        subprocess.run(["git", "config", "diff.sqlite.command", f"{sys.argv[0]} diff"])
        subprocess.run(["git", "config", "merge.sqlite.name", "sqlite merge"])
        subprocess.run(["git", "config", "merge.sqlite.driver", f"{sys.argv[0]} merge %O %A %B %L %P"])

        # git show/apply will not use the external diff program by default
        subprocess.run(["git", "config", "alias.show-sql", "show --ext-diff -m"])
        print("'git show-sql' alias added. Use it instead of 'git show' for sqlite files.")
        print("This script bakes the absolute path to this repo into the repo config. If you move the folder you must rerun this script.")

