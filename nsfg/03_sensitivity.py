exec(open('02_analysis.py').read().split('for sex,g in a.groupby')[0])
for sex,g in a.groupby('sex'):
 print(sex)
 for y in ['current','ever','union']:
  h=g[g.LIFPRTNR>0]; print('positive',y,ols(h,f'{y}~logn+C(AGER)+C(HIEDUC)+C(HISPRACE2)+C(wave)').loc['logn'].round(4).to_dict())
  h=g[g.OPPLIFENUM.between(1,50)].copy();h['logn']=np.log2(1+h.OPPLIFENUM);print('broad positive',y,ols(h,f'{y}~logn+C(AGER)+C(HIEDUC)+C(HISPRACE2)+C(wave)').loc['logn'].round(4).to_dict())
  h=g[g.ever==1];
  if y=='current':print('ever married current',ols(h,f'{y}~logn+C(AGER)+C(HIEDUC)+C(HISPRACE2)+C(wave)').loc['logn'].round(4).to_dict())
# additional variables
exec(open('01_import.py').read().split('cols=')[0])
ls=[]
for wave in ['2015_2017','2017_2019']:
 h=read(next(p.glob(f'{wave}_FemRespSetup*')),next(p.glob(f'{wave}_FemRespData*')),['CASEID','PTSB4MAR','FMAR1AGE','MAREND01','MARDAT01','MARDIS01']);h['wave']=wave;ls.append(h)
f=a[a.sex=='F'].merge(pd.concat(ls),on=['wave','CASEID'],validate='one_to_one');f=f[f.ever==1].copy();print('premarital counts', f.PTSB4MAR.value_counts(dropna=False).sort_index().to_dict());print('end1',f.MAREND01.value_counts(dropna=False).to_dict());print('mar age',f.FMAR1AGE.value_counts(dropna=False).sort_index().to_dict());print('year',f.MARDAT01.min(),f.MARDAT01.max())
f=f[f.PTSB4MAR.between(0,20)].copy();f['premcat']=pd.cut(f.PTSB4MAR,[-1,0,1,4,9,20],labels=['0','1','2-4','5-9','10+']);print('pre N',len(f));
for c,h in f.groupby('premcat',observed=True):print(c,len(h),'mean age at first marriage',np.average(h.FMAR1AGE,weights=h.w),'currently married',np.average(h.current,weights=h.w),'first-end dist',h.MAREND01.value_counts(dropna=False).to_dict())
