"""Weighted prevalence in the 25–49 reference population, from original NSFG files."""
from pathlib import Path
import pandas as pd,numpy as np,re
from scipy.stats import t
p=Path('input')
def read(dct,dat,cols):
 a={m[2].upper():(int(m[1])-1,int(m[3])) for m in re.finditer(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)(?:\.\d+)?f',dct.read_text(errors='replace'))}
 assert set(cols)<=set(a),set(cols)-set(a)
 return pd.read_fwf(dat,colspecs=[(a[c][0],sum(a[c])) for c in cols],names=cols)
frames=[]
for sex,kind,wt in [('F','FemResp','Female'),('M','Male','Male')]:
 w=read(next(p.glob(f'*{wt}WgtSetup*')),next(p.glob(f'*{wt}WgtData*')),['CASEID','WGT2015_2019'])
 for wave in ['2015_2017','2017_2019']:
  d=read(next(p.glob(f'{wave}_{kind}Setup*')),next(p.glob(f'{wave}_{kind}Data*')),['CASEID','AGER','HADSEX','LIFPRTNR','OPPLIFENUM','SEST','SECU']).merge(w,on='CASEID',validate='one_to_one',how='left');d['sex']=sex;frames.append(d)
D=pd.concat(frames,ignore_index=True).rename(columns={'WGT2015_2019':'w'});idx=pd.MultiIndex.from_frame(D[['SEST','SECU']].drop_duplicates().sort_values(['SEST','SECU']));df=len(idx)-D.SEST.nunique();crit=t.ppf(.975,df)
def estimate(g,y):
 w=g.w.to_numpy();x=y.to_numpy(float);b=w@x/w.sum();z=pd.DataFrame({'stratum':g.SEST.to_numpy(),'psu':g.SECU.to_numpy(),'u':w*(x-b)/w.sum()}).groupby(['stratum','psu']).u.sum().reindex(idx,fill_value=0);v=sum(len(s)/(len(s)-1)*((s-s.mean())**2).sum() for _,s in z.groupby(level=0));se=np.sqrt(v);return b,b-crit*se,b+crit*se
rows=[]
for sex,thresholds in [('F',[15,17,23,24]),('M',[15,40,48])]:
 a=D[D.AGER.between(25,49)&(D.sex==sex)&D.HADSEX.isin([1,2])];rows.append(dict(sex=sex,measure='HADSEX',threshold='never',n=len(a),**dict(zip(['estimate','lo','hi'],estimate(a,a.HADSEX==2)))))
 for var in ['OPPLIFENUM','LIFPRTNR']:
  for n in thresholds:
   g=a[a[var].between(1,50)]
   rows.append(dict(sex=sex,measure=var,threshold=n,n=len(g),**dict(zip(['estimate','lo','hi'],estimate(g,g[var]>=n)))))
pd.DataFrame(rows).to_csv('threshold_estimates.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
