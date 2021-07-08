import sqlite3

EFETCH = "esearch -db 'assembly' -query "${QUERY}" | "

conn = sqlite3.connect('GIs.sqlite')
conn.row_factory = sqlite3.Row

seq = []
for row in conn.execute(f'''SELECT gbuid FROM sources as src JOIN strain as st ON src.strain = st.id WHERE seq = NULL AND gbuid;'''):
    seq.append({name: row[""]})