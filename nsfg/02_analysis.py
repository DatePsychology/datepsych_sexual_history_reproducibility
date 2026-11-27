import pandas as pd,numpy as np
from scipy.stats import t
import patsy
D=pd.read_pickle('nsfg.pkl');D['current']=(D.FMARITAL==1).astype(int);D['ever']=(D.FMARNO>0).astype(int);D['union']=D.RMARITAL.isin([1,2]).astype(int);D['former']=D.FMARITAL.isin([2,3,4]).astype(int);D['cat']=pd.cut(D.LIFPRTNR,[-1,0,1,4,9,14,50],labels=['0','1','2-4','5-9','10-14','15+']);D['logn']=np.log2(1+D.LIFPRTNR)
idx=pd.MultiIndex.from_frame(D[['SEST','SECU']].drop_duplicates().sort_values(['SEST','SECU']));df=len(idx)-len(idx.levels[0]); crit=t.ppf(.975,df)
def v(scores,dat):
 a=pd.DataFrame(scores,index=dat.index).assign(SEST=dat.SEST,SECU=dat.SECU).groupby(['SEST','SECU']).sum().reindex(idx,fill_value=0); out=np.zeros((a.shape[1],a.shape[1]));
 for _,x in a.groupby(level=0):
  z=x.to_numpy();z=z-z.mean(axis=0);out+=len(z)/(len(z)-1)*(z.T@z)
 return out
def diff(a,x,y):
 scores=[];means=[];ns=[]
 for i in [0,1]:
  q=a[y]==i;w=a.w*q; mu=(w*a[x]).sum()/w.sum();means.append(mu);ns.append(q.sum());scores.append(w*(a[x]-mu)/w.sum())
 b=means[0]-means[1];se=v((scores[0]-scores[1]).to_numpy()[:,None],a)[0,0]**.5;return ns,means,b,b-crit*se,b+crit*se,2*t.sf(abs(b/se),df)
def ols(a,formula):
 Y,X=patsy.dmatrices(formula,a,return_type='dataframe');a=a.loc[X.index];x=X.to_numpy();y=Y.to_numpy().ravel();w=a.w.to_numpy(); inv=np.linalg.inv(x.T@(w[:,None]*x)); b=inv@(x.T@(w*y));scores=x*(w*(y-x@b))[:,None];V=inv@v(scores,a)@inv; se=np.sqrt(np.diag(V));return pd.DataFrame({'b':b,'lo':b-crit*se,'hi':b+crit*se,'p':2*t.sf(abs(b/se),df)},index=X.columns)
a=D[D.AGER.between(30,49)].copy()
for sex,g in a.groupby('sex'):
 print('\nSEX',sex, len(g))
 for x in ['LIFPRTNR','OPPLIFENUM']:
  q=g[g[x].between(0,50)].copy(); print(x,'current',diff(q,x,'current'),'ever',diff(q,x,'ever'),'union',diff(q,x,'union'))
 for state,h in g.groupby('FMARITAL'):
  print('state',state,len(h),np.average(h.LIFPRTNR,weights=h.w))
 for c,h in g.groupby('cat',observed=True):
  print('cat',c,len(h),*[round(np.average(h[k],weights=h.w)*100,2) for k in ['current','ever','union','former']])
 for y in ['current','ever','union']:
  print(y,'log slope age',ols(g,f'{y}~logn+C(AGER)').loc['logn'].round(4).to_dict(),'adj',ols(g,f'{y}~logn+C(AGER)+C(HIEDUC)+C(HISPRACE2)+C(wave)').loc['logn'].round(4).to_dict())
  print('cat adj',ols(g,f'{y}~C(cat, Treatment(reference="2-4"))+C(AGER)+C(HIEDUC)+C(HISPRACE2)+C(wave)').filter(like='cat',axis=0).round(4).to_string())
 for lo,hi in [(30,34),(35,39),(40,49)]:
  h=g[g.AGER.between(lo,hi)];print('age',lo,hi,'current',diff(h,'LIFPRTNR','current')[2:5],'ever',diff(h,'LIFPRTNR','ever')[2:5])
for y in ['current','ever','union']:
 print('interaction',y,ols(a,f'{y}~logn*C(sex)+C(AGER)*C(sex)+C(HIEDUC)*C(sex)+C(HISPRACE2)*C(sex)+C(wave)*C(sex)').loc['logn:C(sex)[T.M]'].round(4).to_dict())
# Direct comparison in means and test sex difference of current-marital gap using design
s=[]
for sex in ['F','M']:
 for cur in [0,1]:
  q=(a.sex==sex)&(a.current==cur);w=a.w*q; mu=(w*a.LIFPRTNR).sum()/w.sum();s.append(w*(a.LIFPRTNR-mu)/w.sum())
bf=diff(a[a.sex=='F'],'LIFPRTNR','current')[2];bm=diff(a[a.sex=='M'],'LIFPRTNR','current')[2];se=v((s[2]-s[3]-s[0]+s[1]).to_numpy()[:,None],a)[0,0]**.5;print('sex diff gaps',bm-bf,bm-bf-crit*se,bm-bf+crit*se,2*t.sf(abs((bm-bf)/se),df))
