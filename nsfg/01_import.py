import pandas as pd,numpy as np,re,glob
from scipy.stats import t
from pathlib import Path
p=Path('input')
def read(dct,dat,cols):
 a={m[2].upper():(int(m[1])-1,int(m[3])) for m in re.finditer(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)(?:\.\d+)?f',dct.read_text(errors='replace'))}
 assert set(cols)<=set(a),set(cols)-set(a)
 return pd.read_fwf(dat,colspecs=[(a[c][0],sum(a[c])) for c in cols],names=cols)
cols='CASEID AGER FMARITAL RMARITAL FMARNO LIFPRTNR OPPLIFENUM SEST SECU HIEDUC HISPRACE2 RELDLIFE'.split()
all=[]
for sex,kind,wt in [('F','FemResp','Female'),('M','Male','Male')]:
 w=read(next(p.glob(f'*{wt}WgtSetup*')),next(p.glob(f'*{wt}WgtData*')),['CASEID','WGT2015_2019'])
 for wave in ['2015_2017','2017_2019']:
  d=read(next(p.glob(f'{wave}_{kind}Setup*')),next(p.glob(f'{wave}_{kind}Data*')),cols).merge(w,on='CASEID',validate='one_to_one',how='left');d['sex']=sex;d['wave']=wave;all.append(d)
d=pd.concat(all,ignore_index=True).rename(columns={'WGT2015_2019':'w'})
print(d.shape,d.w.isna().sum(),d.groupby(['sex','wave']).size());print('strata',d.groupby('wave').SEST.nunique(),d.groupby(['wave','SEST']).SECU.nunique().value_counts(),d.groupby('SEST').SECU.nunique().value_counts());
for c in cols[2:]:print(c,d[c].value_counts(dropna=False).sort_index().to_dict() if c not in ['SEST','OPPLIFENUM','LIFPRTNR'] else (d[c].min(),d[c].max(),d[c].value_counts(dropna=False).sort_index().tail(20).to_dict()))
d.to_pickle('nsfg.pkl')
