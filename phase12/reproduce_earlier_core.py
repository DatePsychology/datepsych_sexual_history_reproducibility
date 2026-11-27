"""Independent reconstruction of key earlier NLSY estimates; not original Phase 7/8 source.
Retains the Phase 9 landmark construction. Uses a fresh seeded bootstrap, so
percentile intervals need not match the archived bootstrap stream exactly.
"""
from pathlib import Path
import numpy as np,pandas as pd,patsy,statsmodels.api as sm
root=Path(__file__).resolve().parent.parent;target=root/'phase12'
source=(root/'phase9'/'phase9.py').read_text().split("base='age97 + duration",1)[0]
ns={'__file__':str(root/'phase9'/'phase9.py')};exec(compile(source,'phase9_setup','exec'),ns);land=ns['land'];panel=ns['panel'];common=ns['common']
# Conditional fixed-age20 and matched age20/age25 KM point estimates.
def km(g,horizon):
 tt=np.minimum(g.time.to_numpy(int),horizon+1);ee=(g.event.to_numpy()==1)&(tt<=horizon);w=g.weight.to_numpy(float)
 at=np.cumsum(np.bincount(tt,weights=w,minlength=horizon+2)[::-1])[::-1][:horizon+1]
 ev=np.bincount(tt[ee],weights=w[ee],minlength=horizon+2)[:horizon+1];ok=at[1:]>0
 return 1-np.prod(1-ev[1:][ok]/at[1:][ok])
def rd(g,ex,horizon):
 return km(g[g[ex]=='10+'],horizon)-km(g[g[ex]=='1-4'],horizon)
a20=land[land.landmark==20];a25=land[land.landmark==25];m=a25.merge(a20[['id','maximum','cumulative']],on='id',suffixes=('','_age20'),validate='one_to_one');assert len(m)==4620
rng=np.random.default_rng(20261127);rows=[]
for ex in ['maximum','cumulative']:
 for sex,g in m.groupby('sex'):
  old=rd(g,ex+'_age20',180);new=rd(g,ex,180);draws=[]
  for _ in range(1500):
   b=g.iloc[rng.integers(0,len(g),size=len(g))];draws.append(rd(b,ex,180)-rd(b,ex+'_age20',180))
  lo,hi=np.quantile(draws,[.025,.975]);rows.append(dict(exposure=ex,sex=sex,n=len(g),fixed_age20_rd=old,updated_age25_rd=new,attenuation=new-old,lo=lo,hi=hi))
pd.DataFrame(rows).to_csv(target/'paired_attenuation_independent.csv',index=False)
# Replicate the central current and lagged coefficient sequence from Phase 8.
rows=[];audits=[]
for label,d,tag,offset in [('original_style_current',panel.dropna(subset=['current_recent','current_cumulative']),'current',False),('common_current',common,'current',True),('common_lagged',common,'lag',True),('annual_interval_current',common[common.months.between(6,18)],'current',True),('annual_interval_lagged',common[common.months.between(6,18)],'lag',True)]:
 base='age97 + duration + C(race4)';recent=f'C({tag}_recent, Treatment(reference="1"))';hist=f'{tag}_cumulative_log + {tag}_cumulative_zero';expand=' + C(mother_ed) + C(father_ed) + C(income_group) + C(cohab) + C(prior_cohab)'
 for sex,g in d.groupby('sex'):
  for name,f in [('M2_history',base+' + '+hist),('M3_both',base+' + '+recent+' + '+hist),('M6_expanded_background',base+' + '+recent+' + '+hist+expand)]:
   X=patsy.dmatrix(f,g,return_type='dataframe',NA_action='raise');rank=np.linalg.matrix_rank(X);assert rank==X.shape[1];off=np.log(g.months.to_numpy()/12) if offset else np.zeros(len(g));fit=sm.GLM(g.event.to_numpy(),X,family=sm.families.Binomial(sm.families.links.CLogLog()),offset=off).fit(cov_type='cluster',cov_kwds={'groups':g.id.to_numpy()},maxiter=150);assert fit.converged;audits.append(dict(analysis=label,sex=sex,model=name,intervals=len(g),rank=rank,parameters=X.shape[1]));ci=fit.conf_int()
   for term,b in fit.params.items():rows.append(dict(analysis=label,sex=sex,specification=name,term=term,coefficient=b,standard_error=fit.bse[term],lo=ci.loc[term,0],hi=ci.loc[term,1]))
pd.DataFrame(rows).to_csv(target/'event_history_independent.csv',index=False);pd.DataFrame(audits).to_csv(target/'earlier_core_audit.csv',index=False)
print('Independent earlier-core checks completed.')
