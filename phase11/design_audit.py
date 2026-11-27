from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,patsy
from scipy.stats import t as student
root=Path(__file__).resolve().parent.parent;out=Path(__file__).resolve().parent
s=(out/'audit.py').read_text();s=s[:s.index('res=[];background=')];exec(compile(s,str(out/'audit.py'),'exec'),globals());out=root/'phase11'
design=pd.DataFrame(np.load(out/'_design_private.npz')['design'],columns=['id','stratum','psu']);assert len(design)==8984 and design.id.is_unique and set(design.id)==set(data.R0000100);assert set(design.psu)=={1,2};assert design.groupby('stratum').psu.nunique().eq(2).all();idx=pd.MultiIndex.from_frame(design[['stratum','psu']].drop_duplicates().sort_values(['stratum','psu']));df=len(idx)-design.stratum.nunique();print('design strata',design.stratum.nunique(),'PSU',len(idx),'df',df,flush=True)
a=a.merge(design,on='id',validate='one_to_one');crit=student.ppf(.975,df);background=' + C(birth_year)+C(race4)+C(mother_ed)+C(father_ed)+C(income_group)';results=[]
def variance(s):
 z=pd.DataFrame(s).groupby(level=[0,1]).sum().reindex(idx,fill_value=0);v=0
 for _,g in z.groupby(level=0):
  g=g.to_numpy();v+=len(g)/(len(g)-1)*np.sum((g-g.mean(axis=0))**2,axis=0)
 return v
for ex in ['maximum','cumulative','rawmaxgroup','regularmaxgroup']:
 d=a[a[ex].isin(['1-4','10+'])].copy();d['high']=(d[ex]=='10+').astype(int);cuts=[0,24,60,120,180,240];p=[]
 for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
  q=d[d.time>l].copy();q['months']=np.minimum(q.time,r)-l;q['failure']=((q.event==1)&(q.time<=r)).astype(int);q['period']=str(j);p.append(q)
 pp=pd.concat(p,ignore_index=True);X=patsy.dmatrix('0+C(period):C(high):C(sex)'+background,pp,return_type='dataframe',NA_action='raise');x=X.to_numpy();w=pp.weight.to_numpy(float);w/=w.mean();fit=sm.GLM(pp.failure,X,family=sm.families.Poisson(),offset=np.log(pp.months),freq_weights=w).fit(maxiter=100);assert fit.converged and np.linalg.matrix_rank(x)==x.shape[1];mu=fit.fittedvalues.to_numpy();beta=fit.params.to_numpy();bread=np.linalg.pinv(x.T@((w*mu)[:,None]*x));score=pd.DataFrame(x*(w*(pp.failure.to_numpy()-mu))[:,None],index=pp.id).groupby(level=0).sum();inf=score.to_numpy()@bread.T;inf=pd.DataFrame(inf,index=score.index);joint={};points={}
 for sex,g in d.groupby('sex'):
  wg=g.weight.to_numpy()/g.weight.sum();pred=[]
  for high in [0,1]:
   q=g.copy();q['high']=high;H=np.zeros(len(q));dH=np.zeros((len(q),len(beta)))
   for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
    q['period']=str(j);xx=np.asarray(patsy.build_design_matrices([X.design_info],q)[0]);h=(r-l)*np.exp(xx@beta);H+=h;dH+=h[:,None]*xx
   S=np.exp(-H);pred.append((1-wg@S,(wg*S)@dH,1-S))
  low,gl,pl=pred[0];high,gh,ph=pred[1];rd=high-low;grad=gh-gl;model=pd.Series(inf.to_numpy()@grad,index=inf.index);emp=pd.Series(wg*(ph-pl-rd),index=g.id);total=model.add(emp,fill_value=0);designinf=total.to_frame('score').join(design.set_index('id')).set_index(['stratum','psu']).score;se=np.sqrt(variance(designinf)[0]);joint[sex]=designinf.groupby(level=[0,1]).sum().reindex(idx,fill_value=0);points[sex]=rd;results.append(dict(exposure=ex,sex=sex,n=len(g),n_high=int(g.high.sum()),rd=rd,se=se,lo=rd-crit*se,hi=rd+crit*se,df=df))
 diff=points['Women']-points['Men'];se=np.sqrt(variance(joint['Women']-joint['Men'])[0]);results.append(dict(exposure=ex,sex='Women minus Men',n=len(d),rd=diff,se=se,lo=diff-crit*se,hi=diff+crit*se,df=df))
pd.DataFrame(results).to_csv(out/'survey_design_estimates.csv',index=False);print(pd.DataFrame(results).round(4).to_string(index=False))
