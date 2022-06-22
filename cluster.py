#!/usr/bin/env python3
"""
Generate clusters.html
"""
#import os
#import gffutils
#import functools
#import sqlite3

identity_threshold = 90

if __name__ == '__main__':
  clusters = []
  with open('blastn_coverage_filtered_paths.tabular') as f:
    for line in f:
      r = line.split("\t")
      if int(r[6]) < identity_threshold: continue
      for c in clusters:
        if r[0] in c or r[1] in c:
          c[r[0]] = r[4]
          c[r[1]] = r[5]
          break
      else:
        clusters.append({r[0]: r[4], r[1]: r[5]})

  i = 0
  while i < len(clusters):
    a = i + 1
    while a < len(clusters):
      if clusters[a] and set(clusters[i].keys()) & set(clusters[a].keys()):
        clusters[i] |= clusters[a]
        del clusters[a]
      else:
        a += 1
    i += 1

  for i, c in enumerate(clusters):
    print(f"{i}.html", '/mnt/'+' /mnt/'.join(c.values()))

  #conn = sqlite3.connect('GIs.sqlite')
  #conn.row_factory = sqlite3.Row

  #names = {id: name for id, name in conn.execute(f'''SELECT id, name from genomic_islands''')}

  #with open('clusters.html', 'w') as f:
  #  print("""<html><body>
  #  <style>td { border-bottom: 1px solid black;} img { width: 1400px; }</style>
  #  <table>
  #  """, file=f)
  #  for c in clusters:
  #    print("<tr>", file=f)
  #    k = [int(a) for a in c.keys()]
  #    k.sort()
  #    print('<td>' + ' '.join(f'<span title="{names.get(int(a))}">{str(a)}</span>' for a in k) + '</td><td>\n' + '\n'.join(
  #      '<img src="prokka/' + os.path.splitext(os.path.basename(a))[0] + '.svg" />' for a in c.values()) + '</td>', file=f)
  #    cds = {}
  #    for k,a in c.items():
  #      fname = f'prokka/{os.path.splitext(os.path.basename(a))[0]}.gff'
  #      if os.path.exists(fname):
  #        db = gffutils.create_db(fname, ':memory:')
  #        cds[k] = set((feat.attributes.get('gene') or feat.attributes.get('Name') or feat.attributes.get('ID'))[0] for feat in db.features_of_type('CDS'))
  #      else:
  #        cds[k] = None
  #    cdsvalues = list(z for z in cds.values() if z)
  #    commoncds = functools.reduce(lambda x, y: x & y, cdsvalues) if len(cdsvalues) else set()
  #    for k,a in cds.items():
  #      if a and len(a):
  #        print(f"<td>{k}<ul><li>{'</li><li>'.join(a - commoncds)}</li></td>", file=f)
  #      else:
  #        print(f"<td>gff not found for {k}</td>", file=f)
  #    print("</tr>", file=f)
  #  print("""</table>
  #  </body></html>
  #  """, file=f)