import zipfile,numpy as np
from pathlib import Path
root=Path(__file__).resolve().parent.parent
rows=[]
with zipfile.ZipFile(root/'nlsy97_all_1997-2023.zip') as z,z.open('nlsy97_all_1997-2023.csv') as f:
 h=[v.strip(b'"') for v in f.readline().strip().split(b',')];idx=[h.index(v) for v in [b'R0000100',b'R1489700',b'R1489800']]
 for i,line in enumerate(f):
  a=line.split(b',',max(idx)+1);rows.append([int(a[j].strip(b'"')) for j in idx])
  if i%1000==0: print(i,flush=True)
 print('done',i+1,flush=True)
np.savez(root/'phase11'/'_design_private.npz', design=np.array(rows,dtype=int))
