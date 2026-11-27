from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter
R=Path(__file__).resolve().parent.parent;O=Path(__file__).resolve().parent/'figures';O.mkdir(exist_ok=True)
# Regenerate the unchanged figures using the original figure program; its paths are portable in this archive.
s=(Path(__file__).resolve().parent/'make_figures.py').read_text().replace("R=Path('/mnt/data'); O=R/'substack_draft'/'figures'", "R=Path(__file__).resolve().parent.parent; O=Path(__file__).resolve().parent/'figures'")
exec(compile(s,'make_figures.py','exec'),globals())
# Figure 2: target vs reporter labels, with fewer repeated labels.
fig,ax=plt.subplots(figsize=(8.4,4.6));labels=['17+\nmedian','23+\nmean','40+\nmedian','48+\nrounded mean'];vals=[10.27,6.49,8.27,6.78];x=np.arange(4);p=ax.bar(x,vals,color=[F,F,M,M],width=.65);ax.bar_label(p,labels=[f'{v:.1f}%' for v in vals],padding=4);ax.set_xticks(x,labels);ax.set_ylim(0,14);ax.set_ylabel('Weighted population share (%)');ax.grid(axis='y',alpha=.16);ax.set_title('How uncommon are the rejection-case counts?',loc='left',pad=13);ax.text(.5,13.3,'Female partners\nreported by men',ha='center',va='top',fontsize=10,color=F);ax.text(2.5,13.3,'Male partners\nreported by women',ha='center',va='top',fontsize=10,color=M);save(fig,'02_population_thresholds.png')
# Figure 5: actual design-based confidence intervals for risk differences.
forest('05_nsfg_status.png','Marriage associations depend on the outcome',['Currently married · men','Currently married · women','Ever married · men','Ever married · women','Married/cohabiting · men','Married/cohabiting · women'],[-16.23,-24.18,-7.75,-14.38,-8.33,-15.79],[-21.86,-30.16,-12.88,-19.50,-12.77,-21.70],[-10.61,-18.19,-2.62,-9.26,-3.89,-9.87],[M,F,M,F,M,F],'15+ versus 2–4 (percentage points; 95% CI)',(-34,3))
# Figures 6 and 6b: standardized absolute risks and their own design-based intervals.
q=pd.read_csv(R/'phase12'/'reviewer_sensitivity.csv');q=q[q.estimand=='risk']
def risks(label,groups,name,title):
 fig,ax=plt.subplots(figsize=(8.4,4.8));x=np.arange(len(groups));
 for i,(sex,c) in enumerate([('Men',M),('Women',F)]):
  t=q[(q.analysis==label)&(q.sex==sex)].set_index('group').loc[groups];b=t.estimate.to_numpy()*100;lo=t.lo.to_numpy()*100;hi=t.hi.to_numpy()*100;xx=x+(i-.5)*.25
  ax.errorbar(xx,b,yerr=[b-lo,hi-b],fmt='o',capsize=4,color=c,label=sex,lw=2,ms=7)
  for j,v in enumerate(b):ax.annotate(f'{v:.1f}%',(xx[j],v),xytext=(0,10 if i else -15),textcoords='offset points',ha='center',fontsize=9,color=c)
 ax.set_xticks(x,[g.replace('-','–')+'\n'+f"n={int(q[(q.analysis==label)&(q.group==g)&(q.sex=='Men')].n.iloc[0]):,} men /\n{int(q[(q.analysis==label)&(q.group==g)&(q.sex=='Women')].n.iloc[0]):,} women" for g in groups]);ax.set_ylabel('Estimated first marriage by age 40 (%)');ax.set_ylim(0,85);ax.set_title(title,loc='left',pad=15);ax.grid(axis='y',alpha=.16);ax.legend(frameon=False,loc='lower center',ncol=2);save(fig,name)
risks('primary_two_groups',['1-4','10+'],'07_prospective_risks.png','First marriage by age 40: standardized estimates')
risks('all_categories',['0','1-4','5-9','10+'],'10_all_categories.png','All four age-20 maximum-index categories')
# Figure 9: exponentiated coefficient CIs from the actual fitted Phase 8 models.
z=pd.read_csv(R/'event_history_estimates.csv');fig,axes=plt.subplots(2,1,figsize=(8.3,6.8),sharex=True);cats=['0','2','3','4','5+'];xx=np.arange(5)
for ax,sex in zip(axes,['Men','Women']):
 for i,(lab,tag,c) in enumerate([('Contemporaneous','current',N),('Strictly lagged','lagged',M)]):
  t=z[(z.analysis=='common_'+tag)&(z.specification=='M6_expanded_background')&(z.sex==sex)].copy();t['cat']=t.term.str.extract(r'\[T\.(.*?)\]');t=t[t.term.str.contains('_recent')].set_index('cat').loc[cats];b=np.exp(t.coefficient);lo=np.exp(t.ci95_lower);hi=np.exp(t.ci95_upper);ax.errorbar(xx+(i-.5)*.23,b,yerr=[b-lo,hi-b],fmt='o',lw=2,ms=6,capsize=3,color=c,label=lab)
 ax.axhline(1,color='#8392a0',ls='--',lw=1);ax.set_yscale('log');ax.set_ylim(.09,1.35);ax.yaxis.set_major_locator(FixedLocator([.1,.2,.5,1]));ax.yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:g}'));ax.grid(axis='y',alpha=.16);ax.set_title(sex,loc='left',fontsize=13);ax.set_ylabel('Hazard ratio (log scale)');ax.set_xticks(xx,cats);ax.tick_params(axis='x',length=0)
axes[0].legend(loc='upper center',bbox_to_anchor=(.5,1.27),ncol=2,frameon=False);axes[-1].set_xlabel('Recent partners (reference: exactly one)');fig.suptitle('Prospective measurement attenuates associations',x=.1,ha='left',fontsize=15,fontweight='bold');fig.subplots_adjust(top=.85,hspace=.35);save(fig,'09_timing.png')
