import zipfile,numpy as np
from pathlib import Path
root=Path(__file__).resolve().parent.parent
names=['R0000100','R0002100','R0532200','R1200300','R1217500','R1486900','R1489700','R1489800']
rows=[]
with zipfile.ZipFile(root/'nlsy97_all_1997-2023.zip') as z,z.open('nlsy97_all_1997-2023.csv') as f:
 h=[v.strip(b'"') for v in f.readline().strip().split(b',')];idx=[h.index(v.encode()) for v in names]
 for i,line in enumerate(f):
  a=line.split(b',',max(idx)+1);rows.append([int(a[j].strip(b'"')) for j in idx])
  if i%1000==0: print(i,flush=True)
 print('done',i+1,flush=True)
np.savez(root/'phase11'/'_extra_private.npz', data=np.array(rows,dtype=int), names=names)
