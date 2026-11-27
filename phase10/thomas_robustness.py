from pathlib import Path
import zipfile, pandas as pd, numpy as np, statsmodels.formula.api as smf
from scipy.stats import norm
root=Path(__file__).resolve().parent.parent
with zipfile.ZipFile(root/'b4e8u-osfstorage-archive.zip') as z:
    d=pd.read_csv(z.open('Data/Study 3.csv'))
rows=[]; country=[]
for n in [1,12,15]:
    # Compare matched frequency profiles within respondent, not different profile mixes.
    a=d[[f'P4L{i}' for i in range(1,16)]].to_numpy(float)
    b=d[[f'P36L{i}' for i in range(1,16)]].to_numpy(float)
    delta=a-b; valid=np.isfinite(delta).sum(1)
    p=d.loc[valid>=n,['Sex','Country','Ori','Age','SOIR']].copy()
    p['drop']=np.nansum(delta[valid>=n],1)/valid[valid>=n];p['male']=(p.Sex=='Male').astype(int)
    for orientation in ['all','heterosexual']:
        q=p if orientation=='all' else p[p.Ori=='Heterosexual']
        for spec in ['country_FE','country_FE_age_SOIR','equal_country_sex_cells','exclude_other']:
            x=q.copy(); formula='drop~male+C(Country)'
            if spec=='country_FE_age_SOIR':formula+='+Age+SOIR';x=x.dropna(subset=['Age','SOIR'])
            if spec=='exclude_other':x=x[x.Country!='OTH']
            if spec=='equal_country_sex_cells':
                x=x[x.Country!='OTH'].copy();w=1/x.groupby(['Country','Sex']).male.transform('size');f=smf.wls(formula,x,weights=w).fit(cov_type='HC3')
            else:f=smf.ols(formula,x).fit(cov_type='HC3')
            beta,se=f.params['male'],f.bse['male'];lo,hi=beta+np.array([-1,1])*norm.ppf(.975)*se;lo90,hi90=beta+np.array([-1,1])*norm.ppf(.95)*se
            rows.append(dict(min_pairs=n,orientation=orientation,spec=spec,n=len(x),male_minus_female=beta,lo95=lo,hi95=hi,p=f.pvalues['male'],lo90=lo90,hi90=hi90,equiv025=lo90>-.25 and hi90<.25,equiv05=lo90>-.5 and hi90<.5))
        if n==15:
            for c,g in q.groupby('Country'):
                if g.male.nunique()!=2:continue
                f=smf.ols('drop~male',g).fit(cov_type='HC3');beta=f.params['male'];se=f.bse['male'];country.append(dict(orientation=orientation,country=c,n=len(g),female_drop=g.loc[g.male==0,'drop'].mean(),male_drop=g.loc[g.male==1,'drop'].mean(),male_minus_female=beta,lo95=beta-1.96*se,hi95=beta+1.96*se))
pd.DataFrame(rows).to_csv(root/'phase10'/'thomas_robustness.csv',index=False);pd.DataFrame(country).to_csv(root/'phase10'/'thomas_country_estimates.csv',index=False)
print(pd.DataFrame(rows).round(4).to_string(index=False));print(pd.DataFrame(country).round(3).to_string(index=False))
