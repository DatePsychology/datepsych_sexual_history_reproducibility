from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,patsy
from scipy.stats import norm
root=Path(__file__).resolve().parent.parent; out=Path(__file__).resolve().parent
maps=pd.read_csv(root/'resolved_variables(3).csv').fillna('').to_dict('records')
refs='R0000100 R0536300 R0536401 R0536402 R1482600 R1204500 R1302400 R1302500 R1236100 Z9073300 Z9073100'.split(); cols=set(refs)|{v for m in maps for k,v in m.items() if k!='Unnamed: 0' and v}
data=pd.read_csv(root/'respondents.csv',usecols=lambda c:c in cols)
nan=np.nan; valid=lambda x:np.isfinite(x) and x>=0; mx=lambda x,y:max(x,y) if valid(x) else y
ed=lambda x:'missing' if not 0<=x<=20 else 'below_hs' if x<12 else 'hs' if x==12 else 'some_college' if x<16 else 'college_plus'
inc=lambda x:'missing' if not np.isfinite(x) or x in [-1,-2,-3,-4,-5] else 'low' if x<=20000 else 'lower_middle' if x<=40000 else 'upper_middle' if x<=75000 else 'high'
cohab=lambda x:'yes' if x==1 else 'no' if x==2 else 'unknown'
prior=lambda x,d:'no' if x==-4 or x>d else 'yes' if x>0 else 'unknown'
group=lambda x:'0' if x==0 else '1-4' if 1<=x<5 else '5-9' if 5<=x<10 else '10+' if x>=10 else None
rc=lambda x:str(int(x)) if valid(x) and x<5 else '5+' if valid(x) else None
panel=[];land=[]
for row in data.to_dict('records'):
 g=lambda k:row.get(k,nan); M=g('Z9073300');M=np.inf if M==-4 else M
 if not M>0:continue
 B=(g('R0536402')-1980)*12+g('R0536401');hist=[]
 for m in maps:
  d=g(m['interview']);y,mo=g(m['interview_y']),g(m['interview_m'])
  if not d>0:d=(y-1980)*12+mo if y>=1997 and 1<=mo<=12 else nan
  if d>0:hist.append((int(d),m['Unnamed: 0'],m))
 hist.sort(key=lambda h:h[0])
 if len(hist)<2 or M<=hist[0][0] or (np.isfinite(M) and M>hist[-1][0]):continue
 first,last=hist[0][0],hist[-1][0]
 base=dict(id=int(g('R0000100')),sex='Men' if g('R0536300')==1 else 'Women',birth_year=int(g('R0536402')),age97=(first-B)/12,race4=str(int(g('R1482600'))),mother_ed=ed(g('R1302500')),father_ed=ed(g('R1302400')),income_group=inc(g('R1204500')),weight=g('R1236100'))
 c=ma=old=nan; anchored=False;missing=0; hh=[]
 for d,y,m in hist:
  if y>2015:break
  raw,ever,since,once,direct=[g(m[k]) for k in ['recent','ever','since','more_than_once','lifetime']]
  r=raw if valid(raw) else 0 if raw==-4 and (ever==0 or since==0 and y<=2011) else 1 if raw==-4 and once==0 and (since==1 or 2002<=y<=2005) else nan
  a=direct if direct>0 and y<=2005 else 0 if ever==0 else 1 if once==0 and (ever==1 or since==1 or 2002<=y<=2005) else nan
  ro=raw if valid(raw) else 0 if raw==-4 and (ever==0 or since==0 and y<=2009) else 1 if raw==-4 and y<=2001 and once==0 and since==1 else nan
  ao=direct if direct>0 and y<=2001 else 0 if ever==0 else 1 if y<=2001 and once==0 and (ever==1 or since==1) else nan
  if valid(ao):old=ao
  elif valid(ro):old=old+ro if valid(old) else ro
  if ever==1 or since==1:ma=mx(ma,1)
  if valid(r):ma=mx(ma,r)
  if valid(a) and not(valid(ma) and a<ma):c=a;ma=mx(ma,a);anchored=True;missing=0
  elif valid(r):c=c+r if valid(c) else r
  else:missing+=1
  if valid(ma):c=mx(c,ma)
  hh.append(dict(date=d,wave=y,marstat=g(m['marstat']),recent=rc(r),cumulative=c if c!=0 or anchored and missing==0 else nan,maximum=ma if ma!=0 or anchored and missing==0 else nan,oldrec=ro,oldcum=old,landcum=c,landmax=ma))
 for l,r in zip(hh[:-1],hh[1:]):
  s,t=l['date'],r['date']
  if s>=M:break
  if t<=s:continue
  p=dict(**base,start=s,end=t,duration=(s-first)/12,age_start=(s-B)/12,months=t-s,event=int(s<M<=t),cohab=cohab(l['marstat']),prior_cohab=prior(g('Z9073100'),s))
  for tag,h in [('lag',l),('current',r)]:
   for k in ['recent','cumulative','maximum']:p[tag+'_'+k]=h[k]
  panel.append(p)
  if p['event']:break
 for age in [20,25]:
  t=B+age*12
  if M<=t or last<=t:continue
  h=next((h for h in reversed(hh) if h['date']<=t),None)
  if h is None or h['wave']>2009 or t-h['date']>18 or h['marstat'] not in [1,2] or not valid(h['oldrec']) or not valid(h['oldcum']):continue
  land.append(dict(**base,landmark=age,start_month=t,time=min(M,last)-t,event=int(M<=last),cohab=cohab(h['marstat']),prior_cohab=prior(g('Z9073100'),h['date']),recent=h['recent'],maximum=group(h['landmax']),cumulative=group(h['landcum'])))
panel=pd.DataFrame(panel);land=pd.DataFrame(land)
for tag in ['lag','current']:
 for e in ['maximum','cumulative']:
  x=panel[tag+'_'+e];panel[tag+'_'+e+'_log']=np.log2(x.clip(lower=1));panel[tag+'_'+e+'_zero']=(x==0).astype(int)
common=panel.dropna(subset=['lag_recent','current_recent','lag_cumulative','current_cumulative']).copy()
print('panel',len(panel),'common',len(common),'landmarks',land.groupby(['landmark','sex']).size().to_dict(),flush=True)
print(land.groupby(['landmark','sex','maximum']).size(),flush=True)
base='age97 + duration + C(race4)'; back=base+' + C(mother_ed)+C(father_ed)+C(income_group)'; rel=' + C(cohab)+C(prior_cohab)'
rec=lambda tag:f'C({tag}_recent, Treatment(reference="1"))'
def fit(d,f):
 X=patsy.dmatrix(f,d,return_type='dataframe',NA_action='raise');assert np.linalg.matrix_rank(X)==X.shape[1],f
 model=sm.GLM(d.event.to_numpy(),X,family=sm.families.Binomial(sm.families.links.CLogLog()),offset=np.log(d.months.to_numpy()/12));f=model.fit(maxiter=150,cov_type='cluster',cov_kwds={'groups':d.id.to_numpy()});assert f.converged;return f
bridge=common.merge(land[land.landmark==20][['id','start_month','maximum','cumulative']],on='id',validate='many_to_one');bridge=bridge[bridge.start>=bridge.start_month]
res=[]
for e in ['maximum','cumulative']:
 d=bridge[bridge[e].isin(['1-4','10+'])].copy();d['early_high']=(d[e]=='10+').astype(int)
 for sex,g in d.groupby('sex'):
  dyn=f'lag_{e}_log';print('bridge',e,sex,len(g),g.event.sum(),g.id.nunique(),flush=True)
  for lab,extra in [('early_only',''),('recent',' + '+rec('lag')),('updated',' + '+dyn),('both',' + '+rec('lag')+' + '+dyn),('full',' + '+rec('lag')+' + '+dyn+rel)]:
   f=fit(g,back+' + early_high'+extra)
   for t in ['early_high',dyn]:
    if t in f.params:
     b=f.params[t];se=f.bse[t];res.append(dict(exposure=e,sex=sex,model=lab,term=t,n=g.id.nunique(),intervals=len(g),events=int(g.event.sum()),b=b,se=se,hr=np.exp(b),lo=np.exp(b-1.96*se),hi=np.exp(b+1.96*se),p=f.pvalues[t]))
pd.DataFrame(res).to_csv(out/'bridge_corrected.csv',index=False)
def influence(f,ids):
 b=f.params.to_numpy();return pd.DataFrame(f.model.score_obs(b)@np.linalg.pinv(-f.model.hessian(b)).T).groupby(ids).sum()
paired=[]
for sex,g in common.groupby('sex'):
 for label,bb,extra in [('basic',base,''),('expanded',back,rel),('flexible_age','bs(age_start,df=5,degree=3) + C(birth_year)+C(race4)+C(mother_ed)+C(father_ed)+C(income_group)',rel)]:
  fits={tag:fit(g,bb+' + '+rec(tag)+f' + {tag}_cumulative_log + {tag}_cumulative_zero'+extra) for tag in ['current','lag']}
  infl={tag:influence(f,g.id.to_numpy()) for tag,f in fits.items()}
  for k in ['0','2','3','4','5+','cumulative_log']:
   terms={tag:f'{tag}_cumulative_log' if k=='cumulative_log' else rec(tag)+f'[T.{k}]' for tag in fits};i={tag:fits[tag].params.index.get_loc(terms[tag]) for tag in fits}
   d=fits['lag'].params.iloc[i['lag']]-fits['current'].params.iloc[i['current']];u=infl['lag'].iloc[:,i['lag']]-infl['current'].iloc[:,i['current']];se=np.sqrt((u*u).sum()*len(u)/(len(u)-1))
   paired.append(dict(sex=sex,model=label,term=k,current=fits['current'].params.iloc[i['current']],lag=fits['lag'].params.iloc[i['lag']],difference=d,se=se,lo=d-1.96*se,hi=d+1.96*se,p=2*norm.sf(abs(d/se))))
pd.DataFrame(paired).to_csv(out/'paired_timing.csv',index=False)
print('DONE',flush=True)
# Standardized cumulative risks in the original age-20 cohorts.
fixed=[]; aud=[]
for e in ['maximum','cumulative']:
 d=land[(land.landmark==20)&land[e].isin(['1-4','10+'])].copy();d['high']=(d[e]=='10+').astype(int);d['recent']=d.recent.replace({'2':'2-4','3':'2-4','4':'2-4'})
 cuts=[0,24,60,120,180,240]; pieces=[]
 for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
  q=d[d.time>l].copy();q['months']=np.minimum(q.time,r)-l;q['failure']=((q.event==1)&(q.time<=r)).astype(int);q['period']=str(j);pieces.append(q)
 pp=pd.concat(pieces,ignore_index=True)
 for spec,extra in [('unadjusted',''),('background',' + C(birth_year)+C(race4)+C(mother_ed)+C(father_ed)+C(income_group)'),('relationship',' + C(birth_year)+C(race4)+C(mother_ed)+C(father_ed)+C(income_group)+C(cohab)+C(prior_cohab)+C(recent)')]:
  X=patsy.dmatrix('0+C(period):C(high):C(sex)'+extra,pp,return_type='dataframe',NA_action='raise'); w=pp.weight.to_numpy(dtype=float);w/=w.mean();x=X.to_numpy();assert np.linalg.matrix_rank(x)==x.shape[1]
  f=sm.GLM(pp.failure,X,family=sm.families.Poisson(),offset=np.log(pp.months),freq_weights=w).fit(maxiter=100);assert f.converged
  beta=f.params.to_numpy();mu=f.fittedvalues.to_numpy(); bread=np.linalg.pinv(x.T@((w*mu)[:,None]*x)); scores=pd.DataFrame(x*(w*(pp.failure.to_numpy()-mu))[:,None]).groupby(pp.id.to_numpy()).sum().to_numpy(); cov=bread@(scores.T@scores)@bread
  aud.append(dict(exposure=e,specification=spec,n=len(d),parameters=x.shape[1],rank=np.linalg.matrix_rank(x)))
  for sex,g in d.groupby('sex'):
   pred={};wg=g.weight.to_numpy()/g.weight.sum()
   for high in [0,1]:
    q=g.copy();q['high']=high;H=np.zeros(len(q));dH=np.zeros((len(q),len(beta)));oo=[]
    for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
     q['period']=str(j);xx=np.asarray(patsy.build_design_matrices([X.design_info],q)[0]);h=(r-l)*np.exp(xx@beta);H+=h;dH+=h[:,None]*xx;S=np.exp(-H);oo.append((1-wg@S,(wg*S)@dH))
    pred[high]=oo
   for j,r in enumerate(cuts[1:]):
    low,gl=pred[0][j];high,gh=pred[1][j];diff=high-low;grad=gh-gl;se=np.sqrt(max(grad@cov@grad,0));fixed.append(dict(exposure=e,specification=spec,sex=sex,n=len(g),age=20+r/12,low=low,high=high,rd=diff,lo=diff-1.96*se,hi=diff+1.96*se))
pd.DataFrame(fixed).to_csv(out/'fixed_standardized_curves.csv',index=False);pd.DataFrame(aud).to_csv(out/'standardization_audit.csv',index=False)
# Diagnose gross risk-set discrepancy without selecting a threshold to fit the paper.
rows=[]
for name,d in [('all',panel),('current_complete',panel.dropna(subset=['current_recent','current_cumulative'])),('common',common)]:
 for age in [0,18,20,21,25]:
  for sex,g in d[d.age_start>=age].groupby('sex'):
   rows.append(dict(sample=name,age_min=age,sex=sex,n=g.id.nunique(),intervals=len(g),interval_years=g.months.sum()/12,events=int(g.event.sum())))
pd.DataFrame(rows).to_csv(out/'riskset_audit.csv',index=False)
print('STANDARDIZATION DONE',flush=True)
