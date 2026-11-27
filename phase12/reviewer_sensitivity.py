from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,patsy
from scipy.stats import t as student
root=Path(__file__).resolve().parent.parent;target=root/'phase12';target.mkdir(exist_ok=True)
s=(root/'phase11'/'expanded_audit.py').read_text().split("for ex in ['maximum'",1)[0];ns={'__file__':str(root/'phase11'/'expanded_audit.py')};exec(compile(s,'expanded_setup','exec'),ns)
a=ns['a'];design=ns['design'];idx=ns['idx'];df=ns['df'];background=ns['background'];variance=ns['variance'];crit=student.ppf(.975,df)
counts={}
for row in ns['data'].to_dict('records'):
 g=lambda k:row.get(k,np.nan);birthday=(g('R0536402')-1980)*12+g('R0536401')+240;waves=set()
 for m in ns['maps']:
  y=m['Unnamed: 0'];d=g(m['interview']);yy,mm=g(m['interview_y']),g(m['interview_m'])
  if not d>0:d=(yy-1980)*12+mm if yy>=1997 and 1<=mm<=12 else np.nan
  if d>0 and d<=birthday and y<=2009 and (g(m['recent'])>=10 or (y<=2005 and g(m['lifetime'])>=10)):waves.add(y)
 counts[int(g('R0000100'))]=len(waves)
a['high_waves']=a.id.map(counts);a.groupby(['sex','maximum','high_waves']).size().rename('n').reset_index().to_csv(target/'corroboration_counts.csv',index=False)
rows=[];audits=[]
for label,d in [('all_categories',a[a.maximum.isin(['0','1-4','5-9','10+'])].copy()),('primary_two_groups',a[a.maximum.isin(['1-4','10+'])].copy()),('corroborated_high',a[(a.maximum=='1-4')|((a.maximum=='10+')&(a.high_waves>=2))].copy())]:
 d['grp']=d.maximum;cuts=[0,24,60,120,180,240];p=[]
 for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
  q=d[d.time>l].copy();q['months']=np.minimum(q.time,r)-l;q['failure']=((q.event==1)&(q.time<=r)).astype(int);q['period']=str(j);p.append(q)
 pp=pd.concat(p,ignore_index=True);X=patsy.dmatrix('0+C(period):C(grp):C(sex)'+background,pp,return_type='dataframe',NA_action='raise');x=X.to_numpy();w=pp.weight.to_numpy(float);w/=w.mean();fit=sm.GLM(pp.failure,X,family=sm.families.Poisson(),offset=np.log(pp.months),freq_weights=w).fit(maxiter=150);rank=np.linalg.matrix_rank(x);assert fit.converged and rank==x.shape[1];mu=fit.fittedvalues.to_numpy();beta=fit.params.to_numpy();bread=np.linalg.pinv(x.T@((w*mu)[:,None]*x));score=pd.DataFrame(x*(w*(pp.failure.to_numpy()-mu))[:,None],index=pp.id).groupby(level=0).sum();inf=pd.DataFrame(score.to_numpy()@bread.T,index=score.index);audits.append(dict(analysis=label,n=len(d),rank=rank,parameters=x.shape[1],max_abs_coefficient=abs(beta).max()))
 allsex={}
 for sex,g in d.groupby('sex'):
  wg=g.weight.to_numpy()/g.weight.sum();pred={}
  for grp in sorted(g.grp.unique()):
   q=g.copy();q['grp']=grp;H=np.zeros(len(q));dH=np.zeros((len(q),len(beta)))
   for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
    q['period']=str(j);xx=np.asarray(patsy.build_design_matrices([X.design_info],q)[0]);h=(r-l)*np.exp(xx@beta);H+=h;dH+=h[:,None]*xx
   S=np.exp(-H);point=1-wg@S;grad=(wg*S)@dH;model=pd.Series(inf.to_numpy()@grad,index=inf.index);emp=pd.Series(wg*(1-S-point),index=g.id);total=model.add(emp,fill_value=0);z=total.to_frame('score').join(design.set_index('id')).set_index(['stratum','psu']).score;u=z.groupby(level=[0,1]).sum().reindex(idx,fill_value=0);pred[grp]=(point,u);se=np.sqrt(variance(u)[0]);rows.append(dict(analysis=label,sex=sex,group=grp,estimand='risk',n=int((g.grp==grp).sum()),events=int(((g.grp==grp)&(g.event==1)&(g.time<=240)).sum()),estimate=point,se=se,lo=point-crit*se,hi=point+crit*se))
  allsex[sex]={}
  for grp,(point,u) in pred.items():
   if grp=='1-4':continue
   rd=point-pred['1-4'][0];du=u-pred['1-4'][1];se=np.sqrt(variance(du)[0]);allsex[sex][grp]=(rd,du);rows.append(dict(analysis=label,sex=sex,group=grp,estimand='difference',n=int((g.grp==grp).sum()),estimate=rd,se=se,lo=rd-crit*se,hi=rd+crit*se))
 for grp in allsex['Women']:
  rd=allsex['Women'][grp][0]-allsex['Men'][grp][0];se=np.sqrt(variance(allsex['Women'][grp][1]-allsex['Men'][grp][1])[0]);rows.append(dict(analysis=label,sex='Women minus Men',group=grp,estimand='difference_of_differences',estimate=rd,se=se,lo=rd-crit*se,hi=rd+crit*se))
 print(label,flush=True)
pd.DataFrame(rows).to_csv(target/'reviewer_sensitivity.csv',index=False);pd.DataFrame(audits).to_csv(target/'model_audit.csv',index=False);print(pd.DataFrame(rows).round(4).to_string(index=False))
