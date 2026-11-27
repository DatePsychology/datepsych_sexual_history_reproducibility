"""Published aggregate counts only, conditional on reconstructed denominators."""
from scipy.stats import fisher_exact,t,norm
import numpy as np,pandas as pd
from pathlib import Path
out=Path(__file__).resolve().parent
rows=[]
for label,a,n,b,m in [('Virgin vs nonvirgin',17,61,135,401),('Partnered vs single',80,215,73,247),('Married vs unmarried',38,93,115,369),('Active vs inactive in 3 years',129,380,24,82)]:
 rows.append(dict(comparison=label,n1=n,yes1=a,n2=m,yes2=b,share1=a/n,share2=b/m,p=fisher_exact([[a,n-a],[b,m-b]])[1]))
pd.DataFrame(rows).to_csv(out/'fisher_tests.csv',index=False)
rows=[]
for sex,delta,stat,df in [('Women',3.33,.912,60.27),('Men',5.58,2.4409,232.92)]:
 se=delta/stat;rows.append(dict(sex=sex,difference=delta,se=se,df=df,lo=delta-t.ppf(.975,df)*se,hi=delta+t.ppf(.975,df)*se))
a,b=rows;v1=a['se']**2;v2=b['se']**2;se=np.sqrt(v1+v2);df=(v1+v2)**2/(v1**2/a['df']+v2**2/b['df']);d=b['difference']-a['difference'];rows.append(dict(sex='Men minus women',difference=d,se=se,df=df,lo=d-t.ppf(.975,df)*se,hi=d+t.ppf(.975,df)*se,p=2*t.sf(abs(d/se),df)))
pd.DataFrame(rows).to_csv(out/'welch_reconstruction.csv',index=False)
pd.DataFrame([dict(reference=label,virgin_share=p,required_ratio=(1-p)/p) for label,p in [('NSFG 25-49',.0364),('Reconstructed sample',61/462),('Reconstructed men with ideal',44/310)]]).to_csv(out/'base_rates.csv',index=False)
rows=[]
for label,n1,m1,n2,m2 in [('Relationship',215,4.57,247,4.73),('Marriage',93,4.35,369,4.73),('Sex within 3 years',380,4.59,82,4.95),('Virginity',61,4.39,401,3.67)]:
 rows.append(dict(partition=label,n1=n1,mean1=m1,n2=n2,mean2=m2,implied_mean=(n1*m1+n2*m2)/(n1+n2)))
pd.DataFrame(rows).to_csv(out/'importance_means.csv',index=False)
print('Minimum combined mean with 50 additional scores of 1:',(rows[0]['implied_mean']*462+50)/512)
print('Illustrative probability ever rejected:',1-.98**7)
