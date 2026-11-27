from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,patsy
from scipy.stats import t as student
root=Path(__file__).resolve().parent.parent;out=Path(__file__).resolve().parent
s=(out/'audit.py').read_text();s=s[:s.index('res=[];background=')];exec(compile(s,str(out/'audit.py'),'exec'),globals());out=root/'phase11'
design=pd.DataFrame(np.load(out/'_design_private.npz')['design'],columns=['id','stratum','psu']);assert len(design)==8984 and design.id.is_unique and set(design.id)==set(data.R0000100);assert set(design.psu)=={1,2};assert design.groupby('stratum').psu.nunique().eq(2).all();idx=pd.MultiIndex.from_frame(design[['stratum','psu']].drop_duplicates().sort_values(['stratum','psu']));df=len(idx)-design.stratum.nunique();print('design strata',design.stratum.nunique(),'PSU',len(idx),'df',df,flush=True)
a=a.merge(design,on='id',validate='one_to_one')
extra=np.load(out/'_extra_private.npz');ee=pd.DataFrame(extra['data'],columns=extra['names']);assert len(ee)==8984 and ee.R0000100.is_unique and set(ee.R0000100)==set(data.R0000100);assert np.array_equal(ee.set_index('R0000100').loc[design.id,['R1489700','R1489800']].to_numpy(),design[['stratum','psu']].to_numpy());print('extra counts',ee.R0532200.value_counts().to_dict(),ee.R1200300.value_counts().to_dict(),ee.R1217500.value_counts().to_dict(),ee.R0002100.value_counts().to_dict(),flush=True)
assert ee.R0532200.value_counts().to_dict()=={0:4589,1:4395};assert ee.R1200300.value_counts().to_dict()=={3:3359,2:2050,4:1990,1:1585};assert ee.R1217500.value_counts().to_dict()=={1:6574,0:2030,2:380}
assert ee.R0002100.value_counts().get(1)==2522 and ee.R0002100.value_counts().get(2)==1991 and ee.R0002100.value_counts().get(-2)==113
assert ee.R1486900.between(0,600).sum()==5961
# Keep religious categories broad and transparent; not a Steensland tradition recode.
def religion(x):
 if x<0:return 'missing'
 if x in [25,26,27]:return 'none'
 if x==1:return 'catholic'
 if x==2:return 'baptist'
 if 3<=x<=13:return 'other_protestant'
 return 'other'
ee['religion']=ee.R0002100.map(religion);ee['parent_relig']=pd.cut(ee.R1486900.where(ee.R1486900.between(0,600)),[-1,200,400,599,600],labels=['0-2','2-4','4-<6','6']).astype('object').fillna('missing');ee['bothbio']=ee.R0532200.astype(str);ee['region']=ee.R1200300.astype(str);ee['urban']=ee.R1217500.astype(str);a=a.merge(ee.rename(columns={'R0000100':'id'})[['id','religion','parent_relig','bothbio','region','urban']],on='id',validate='one_to_one');crit=student.ppf(.975,df);background=' + C(birth_year)+C(race4)+C(mother_ed)+C(father_ed)+C(income_group)+C(religion)+C(parent_relig)+C(bothbio)+C(region)+C(urban)';results=[]
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
pd.DataFrame(results).to_csv(out/'expanded_design_estimates.csv',index=False);print(pd.DataFrame(results).round(4).to_string(index=False))
