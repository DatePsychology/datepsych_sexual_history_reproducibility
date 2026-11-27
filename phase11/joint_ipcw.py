from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,patsy
root=Path(__file__).resolve().parent.parent;out=Path(__file__).resolve().parent
s=(out/'expanded_audit.py').read_text();s=s[:s.index("for ex in ['maximum','cumulative','rawmaxgroup','regularmaxgroup']:")];exec(compile(s,str(out/'expanded_audit.py'),'exec'),globals());out=root/'phase11'
results=[];aud=[]
for ex,horizon,older in [('maximum',240,False),('cumulative',240,False),('maximum',216,False),('maximum',240,True)]:
 d=a[a[ex].isin(['1-4','10+'])].copy();d=d[d.birth_year<=1982].copy() if older else d;d=d.reset_index(drop=True);d['high']=(d[ex]=='10+').astype(int);d['ix']=np.arange(len(d));cuts=[0,24,60,120,180,horizon];pieces=[]
 for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
  q=d[d.time>l].copy();q['months']=np.minimum(q.time,r)-l;q['failure']=((q.event==1)&(q.time<=r)).astype(int);q['censored']=((q.event==0)&(q.time<=r)&(q.time<horizon)).astype(int);q['period']=str(j);pieces.append(q)
 pp=pd.concat(pieces,ignore_index=True);x_c=patsy.dmatrix('C(period)+C(high)*C(sex)'+background+' + C(cohab)+C(prior_cohab)',pp,return_type='dataframe',NA_action='raise');w=pp.weight.to_numpy(float);w/=w.mean();assert np.linalg.matrix_rank(x_c)==x_c.shape[1];fc=sm.GLM(pp.censored,x_c,family=sm.families.Poisson(),offset=np.log(pp.months),freq_weights=w).fit(maxiter=100);assert fc.converged;assert np.isfinite(fc.params).all()
 # G just before each observed month. Within a piece constant censoring rate.
 H=np.zeros((len(d),horizon));months=np.arange(1,horizon+1)
 for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
  q=d.copy();q['period']=str(j);xx=np.asarray(patsy.build_design_matrices([x_c.design_info],q)[0]);rate=np.exp(xx@fc.params.to_numpy());H+=rate[:,None]*np.clip(months-1-l,0,r-l)[None,:]
 raw=np.exp(np.minimum(H,30));obs=d.time.to_numpy()[:,None]>=months[None,:];aud.append(dict(exposure=ex,horizon=horizon,older=older,n=len(d),censor_rank=np.linalg.matrix_rank(x_c),censor_parameters=x_c.shape[1],censored=pp.censored.sum(),multiplier_p99=np.quantile(raw[obs],.99),multiplier_max=np.max(raw[obs]),max_abs_coef=np.abs(fc.params).max()))
 X=patsy.dmatrix('0+C(period):C(high):C(sex)'+background,pp,return_type='dataframe',NA_action='raise');assert np.linalg.matrix_rank(X)==X.shape[1]
 for cap in [1,5,10,20]:
  # Cap 1 = ordinary background standardization.
  weights=np.minimum(raw,cap);E=np.empty(len(pp));D=np.empty(len(pp))
  for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
   mask=pp.period==str(j);q=pp.loc[mask];i=q.ix.to_numpy();W=weights[i,l:r];m=months[l:r];E[mask]=np.sum(W*(m[None,:]<=q.time.to_numpy()[:,None]),axis=1);failure=q.failure.to_numpy(bool);D[mask]=np.where(failure,weights[i,np.clip(q.time.to_numpy(int)-1,0,horizon-1)],0)
  assert (E>0).all();f=sm.GLM(D/E,X,family=sm.families.Poisson(),freq_weights=w*E).fit(maxiter=100);assert f.converged and np.isfinite(f.params).all();beta=f.params.to_numpy()
  for sex,g in d.groupby('sex'):
   wg=g.weight.to_numpy()/g.weight.sum();pred=[]
   for high in [0,1]:
    q=g.copy();q['high']=high;ch=np.zeros(len(q))
    for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
     q['period']=str(j);xx=np.asarray(patsy.build_design_matrices([X.design_info],q)[0]);ch+=(r-l)*np.exp(xx@beta)
    pred.append(1-wg@np.exp(-ch))
   results.append(dict(exposure=ex,horizon=horizon,older=older,cap=cap,sex=sex,n=len(g),low=pred[0],high=pred[1],rd=pred[1]-pred[0]))
 print('done',ex,horizon,older,flush=True)
pd.DataFrame(results).to_csv(out/'joint_ipcw_estimates.csv',index=False);pd.DataFrame(aud).to_csv(out/'joint_ipcw_audit.csv',index=False);print(pd.DataFrame(results).round(4).to_string(index=False));print(pd.DataFrame(aud).round(4).to_string(index=False))
