# Run from a directory containing input/ with CDC.zip extracted.
import pandas as pd,numpy as np,re,glob
from scipy.stats import t
from pathlib import Path
p=Path('input')
def read(dct,dat,cols):
 a={m[2].upper():(int(m[1])-1,int(m[3])) for m in re.finditer(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)(?:\.\d+)?f',dct.read_text(errors='replace'))}
 assert set(cols)<=set(a),set(cols)-set(a)
 return pd.read_fwf(dat,colspecs=[(a[c][0],sum(a[c])) for c in cols],names=cols)

import patsy
from scipy.stats import f,t
w=read(next(p.glob('*FemaleWgtSetup*')),next(p.glob('*FemaleWgtData*')),['CASEID','WGT2015_2019'])
ls=[]
for wave in ['2015_2017','2017_2019']:
 d=read(next(p.glob(f'{wave}_FemRespSetup*')),next(p.glob(f'{wave}_FemRespData*')),['CASEID','AGER','MENARCHE','PTSB4MAR','FMARNO','LIFPRTNR','HIEDUC','HISPRACE2','SEST','SECU']); d['wave']=wave;ls.append(d)
D=pd.concat(ls,ignore_index=True).merge(w,on='CASEID',validate='one_to_one').rename(columns={'WGT2015_2019':'w'});D['ever']=(D.FMARNO>0).astype(int);D['pn']=np.where(D.ever==1,D.PTSB4MAR,D.LIFPRTNR);D.loc[~D.pn.between(0,50),'pn']=np.nan;D['six']=(D.pn>=6).astype(int);D['men']=pd.cut(D.MENARCHE,[5,11,12,13,26],labels=['<=11','12','13','14+']);
idx=pd.MultiIndex.from_frame(D[['SEST','SECU']].drop_duplicates().sort_values(['SEST','SECU']));df=len(idx)-len(idx.levels[0]);crit=t.ppf(.975,df)
def v(scores,dat):
 a=pd.DataFrame(scores,index=dat.index).assign(SEST=dat.SEST,SECU=dat.SECU).groupby(['SEST','SECU']).sum().reindex(idx,fill_value=0); out=np.zeros((a.shape[1],a.shape[1]));
 for _,x in a.groupby(level=0):
  z=x.to_numpy();z=z-z.mean(axis=0);out+=len(z)/(len(z)-1)*(z.T@z)
 return out

def fit(a,formula):
 Y,X=patsy.dmatrices(formula,a,return_type='dataframe');a=a.loc[X.index];x=X.to_numpy();y=Y.to_numpy().ravel();w=a.w.to_numpy();inv=np.linalg.inv(x.T@(w[:,None]*x));b=inv@(x.T@(w*y));V=inv@v(x*(w*(y-x@b))[:,None],a)@inv;se=np.sqrt(np.diag(V)); ii=[i for i,k in enumerate(X.columns) if k.startswith('C(men')];z=b[ii];vv=V[np.ix_(ii,ii)];F=z@np.linalg.solve(vv,z)/len(ii);return pd.DataFrame({'b':b,'lo':b-crit*se,'hi':b+crit*se},index=X.columns).iloc[ii],F,f.sf(F,len(ii),df)
for lo,hi in [(15,44),(25,44)]:
 a=D[D.AGER.between(lo,hi)&D.men.notna()&D.pn.notna()].copy(); print('N',lo,hi,len(a),'df',df);print(a.groupby('men',observed=True).apply(lambda h:pd.Series({'n':len(h),'six':np.average(h.six,weights=h.w),'ever':np.average(h.ever,weights=h.w)}),include_groups=False));
 for y in ['six','ever']:
  for controls in ['C(AGER)+C(HISPRACE2)+C(wave)','C(AGER)+C(HISPRACE2)+C(wave)+C(HIEDUC)']:
   print(y,controls,fit(a,f'{y}~C(men, Treatment(reference="12"))+{controls}'))
