from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm,patsy
root=Path(__file__).resolve().parent.parent;out=Path(__file__).resolve().parent
s=(root/'phase9'/'phase9.py').read_text();s=s[:s.index("base='age97 + duration")];exec(compile(s,str(root/'phase9'/'phase9.py'),'exec'),globals());out=root/'phase11'
raw=[]
for row in data.to_dict('records'):
 g=lambda k:row.get(k,np.nan);B=(g('R0536402')-1980)*12+g('R0536401');t=B+240;hs=[]
 for m in maps:
  y=m['Unnamed: 0'];d=g(m['interview']);yy,mm=g(m['interview_y']),g(m['interview_m'])
  if not d>0:d=(yy-1980)*12+mm if yy>=1997 and 1<=mm<=12 else np.nan
  if d>0 and d<=t and y<=2009:hs.append((d,y,m))
 hs.sort(key=lambda x:x[0]);values=[];recent=[];direct=[];pos=[];reg=[];prev=None
 for d,y,m in hs:
  r,a,ever,since,once=[g(m[k]) for k in ['recent','lifetime','ever','since','more_than_once']]
  if valid(r):
   values.append(r);recent.append(r)
   if y==1997 or prev is not None and 9<=d-prev<=15:reg.append(r)
  if a>0 and y<=2005:values.append(a);direct.append(a)
  if r>0 or a>0 or ever==1 or since==1:pos.append((d-B)/12)
  prev=d
 raw.append(dict(id=int(g('R0000100')),rawmax=max(values,default=np.nan),directmax=max(direct,default=np.nan),regularmax=max(reg+direct,default=np.nan),recentmax=max(recent,default=np.nan),ncounts=len(values),first_positive_age=min(pos,default=np.nan),lag=t-hs[-1][0] if hs else np.nan))
a=land[land.landmark==20].merge(pd.DataFrame(raw),on='id',validate='one_to_one');a['recentcat']=a.recent.replace({'2':'2-4','3':'2-4','4':'2-4'});a['firstposcat']=pd.cut(a.first_positive_age,[-np.inf,15,17,19,np.inf],labels=['<=15','15-17','17-19','>19']).astype(str);a['lagcat']=pd.cut(a.lag,[-1,6,12,18],labels=['0-6','7-12','13-18']).astype(str);a['nc_cat']=pd.cut(a.ncounts,[-1,2,4,6,100],labels=['0-2','3-4','5-6','7+']).astype(str)
a.groupby(['sex','maximum','recentcat'],dropna=False).size().rename('n').reset_index().to_csv(out/'support_audit.csv',index=False)
for c in ['rawmax','regularmax','recentmax','directmax']:a[c+'group']=a[c].map(group)
a.groupby(['sex','maximum','rawmaxgroup'],dropna=False).size().rename('n').reset_index().to_csv(out/'raw_reclassification.csv',index=False)
res=[];background=' + C(birth_year)+C(race4)+C(mother_ed)+C(father_ed)+C(income_group)'
def fitrisk(d,ex,label,extra):
 d=d[d[ex].isin(['1-4','10+'])].copy();d['high']=(d[ex]=='10+').astype(int);cuts=[0,24,60,120,180,240];p=[]
 for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
  q=d[d.time>l].copy();q['months']=np.minimum(q.time,r)-l;q['failure']=((q.event==1)&(q.time<=r)).astype(int);q['period']=str(j);p.append(q)
 pp=pd.concat(p,ignore_index=True);X=patsy.dmatrix('0+C(period):C(high):C(sex)'+extra,pp,return_type='dataframe',NA_action='raise');x=X.to_numpy();rank=np.linalg.matrix_rank(x)
 if rank!=x.shape[1]:print('rank fail',label,rank,x.shape[1]);return
 w=pp.weight.to_numpy(float);w/=w.mean();f=sm.GLM(pp.failure,X,family=sm.families.Poisson(),offset=np.log(pp.months),freq_weights=w).fit(maxiter=100);mu=f.fittedvalues.to_numpy();beta=f.params.to_numpy();bread=np.linalg.pinv(x.T@((w*mu)[:,None]*x));scores=pd.DataFrame(x*(w*(pp.failure.to_numpy()-mu))[:,None]).groupby(pp.id.to_numpy()).sum().to_numpy();cov=bread@(scores.T@scores)@bread
 for sex,g in d.groupby('sex'):
  wg=g.weight.to_numpy()/g.weight.sum();pred=[]
  for high in [0,1]:
   q=g.copy();q['high']=high;H=np.zeros(len(q));dH=np.zeros((len(q),len(beta)))
   for j,(l,r) in enumerate(zip(cuts[:-1],cuts[1:])):
    q['period']=str(j);xx=np.asarray(patsy.build_design_matrices([X.design_info],q)[0]);h=(r-l)*np.exp(xx@beta);H+=h;dH+=h[:,None]*xx
   S=np.exp(-H);pred.append((1-wg@S,(wg*S)@dH))
  low,gl=pred[0];high,gh=pred[1];rd=high-low;grad=gh-gl;se=np.sqrt(max(grad@cov@grad,0));res.append(dict(label=label,exposure=ex,sex=sex,n_low=sum(g.high==0),n_high=sum(g.high==1),low=low,high=high,rd=rd,lo=rd-1.96*se,hi=rd+1.96*se,rank=rank,parameters=x.shape[1],converged=f.converged))
 print(label,flush=True)
fitrisk(a,'maximum','original_background',background)
fitrisk(a,'maximum','cohab_only',background+' + C(cohab)+C(prior_cohab)')
fitrisk(a[a.recentcat!='5+'],'maximum','common_recent_support',background+' + C(cohab)+C(prior_cohab)+C(recentcat)')
fitrisk(a,'rawmaxgroup','direct_report_only',background)
fitrisk(a,'regularmaxgroup','regular_interval_or_lifetime',background)
fitrisk(a,'recentmaxgroup','interval_reports_only',background)
fitrisk(a,'directmaxgroup','lifetime_reports_only_selected',background)
fitrisk(a,'maximum','report_opportunity_adjustment',background+' + C(nc_cat)+C(lagcat)+C(firstposcat)')
fitrisk(a[a.cohab=='no'],'maximum','noncohabiting_background',background+' + C(prior_cohab)')
fitrisk(a[a.birth_year<=1982],'maximum','older_cohort_background',background)
pd.DataFrame(res).to_csv(out/'adversarial_estimates.csv',index=False);print(pd.DataFrame(res).round(4).to_string(index=False))
