"""Extract only the columns needed by the included NLSY scripts, from official ZIP.
All output includes individual records. Keep local; never commit it to GitHub.
"""
from pathlib import Path
import csv,zipfile
import numpy as np,pandas as pd
root=Path(__file__).resolve().parent
maps=pd.read_csv(root/'resolved_variables(3).csv').fillna('')
refs=set(maps.drop(columns='Unnamed: 0').to_numpy().ravel())- {''}
refs.update('R0000100 R0536300 R0536401 R0536402 R1482600 R1204500 R1302400 R1302500 R1236100 Z9073300 Z9073100'.split())
extra='R0000100 R0002100 R0532200 R1200300 R1217500 R1486900 R1489700 R1489800'.split()
rows=[];all_ids=[]
with zipfile.ZipFile(root/'nlsy97_all_1997-2023.zip') as z,z.open('nlsy97_all_1997-2023.csv') as f, (root/'respondents.csv').open('w',newline='') as out:
 header=[v.strip(b'"').decode('ascii') for v in f.readline().strip().split(b',')]
 assert refs.union(extra)<=set(header), 'Required variable missing from source'
 cols=[v for v in header if v in refs];ix=[header.index(v) for v in cols];ei=[header.index(v) for v in extra];stop=max(ix+ei)+1
 writer=csv.writer(out);writer.writerow(cols)
 for i,line in enumerate(f):
  fields=line.split(b',',stop)
  vals=[int(fields[j].strip(b'"\r\n ')) for j in ix];writer.writerow(vals)
  e=[int(fields[j].strip(b'"\r\n ')) for j in ei];rows.append(e);all_ids.append(e[0])
  if (i+1)%1000==0:print(i+1,flush=True)
assert len(all_ids)==len(set(all_ids))==8984
arr=np.array(rows,dtype=int);(root/'phase11').mkdir(exist_ok=True)
np.savez(root/'phase11'/'_extra_private.npz',data=arr,names=extra)
np.savez(root/'phase11'/'_design_private.npz',design=arr[:,[0,6,7]])
print('Extracted 8984 unique IDs; private local files created.')
