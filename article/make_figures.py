from pathlib import Path
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
R=Path(__file__).resolve().parent.parent; O=Path(__file__).resolve().parent/'figures'; O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':15,'axes.titleweight':'bold','axes.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'axes.spines.bottom':False,'xtick.color':'#334155','ytick.color':'#334155','text.color':'#172b3a','axes.labelcolor':'#334155','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','axes.axisbelow':True})
M='#2e6389'; F='#ae6248'; N='#718397'; D='#223a4d'
def save(fig,name):
 fig.savefig(O/name,dpi=180,bbox_inches='tight',pad_inches=.28);plt.close(fig)
def forest(name,title,labels,b,lo,hi,colors,xlabel,limits=None):
 fig,ax=plt.subplots(figsize=(8.5,max(2.9,.51*len(labels)+1.6)));y=np.arange(len(labels))[::-1]; ax.axvline(0,color='#8392a0',ls='--',lw=1);ax.grid(axis='x',alpha=.16)
 for k in range(len(y)):
  ax.errorbar(b[k],y[k],xerr=[[b[k]-lo[k]],[hi[k]-b[k]]],fmt='o',color=colors[k],capsize=3,lw=2,ms=7)
 ax.set_yticks(y,labels);ax.tick_params(axis='y',length=0);ax.set_xlabel(xlabel);ax.set_title(title,loc='left',pad=16);ax.set_ylim(-.6,len(y)-.4)
 if limits: ax.set_xlim(*limits)
 save(fig,name)
def bars(name,title,labels,series,xlabel,ylim=None,percent=False):
 fig,ax=plt.subplots(figsize=(8.2,4.25));x=np.arange(len(labels));w=.72/len(series)
 for j,(label,v,color) in enumerate(series):
  p=ax.bar(x+(j-(len(series)-1)/2)*w,v,w*.93,label=label,color=color,zorder=3)
  ax.bar_label(p,labels=[f'{z:.1f}'+('%' if percent else '') for z in v],padding=4,fontsize=10)
 ax.set_xticks(x,labels);ax.tick_params(axis='x',length=0);ax.grid(axis='y',alpha=.16);ax.set_ylabel(xlabel);ax.set_title(title,loc='left',pad=17)
 if ylim:ax.set_ylim(*ylim)
 if len(series)>1:ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2)
 save(fig,name)
bars('01_denominators.png','Identical subgroup counts, opposite conclusions',['Left','Center','Right'],[('Hypothetical sample A',[104/300*100,27/80*100,22/82*100],M),('Hypothetical sample B',[104/350*100,27/60*100,22/52*100],F)],'Within-group share with no ideal (%)',(0,55),True)
bars('02_population_thresholds.png','How uncommon are the reported rejection-case counts?',['Women: 17+\n(median)','Women: 24+\n(above mean 23*)','Men: 40+\n(median)','Men: 48+\n(above mean 47.5)'],[('Weighted population share',[10.27,6.24,8.27,6.78],M)],'Share at or above threshold (%)',(0,13),True)
# Use 23+ (not >23) for mean threshold. Replace previous chart with unambiguous threshold labels.
bars('02_population_thresholds.png','How uncommon are the reported rejection-case counts?',['Women: 17+\n(median)','Women: 23+\n(mean)','Men: 40+\n(median)','Men: 48+\n(rounded mean)'],[('Weighted population share',[10.27,6.49,8.27,6.78],M)],'Share at or above threshold (%)',(0,13),True)
forest('03_experimental_sex_difference.png','Sex difference in decline in willingness ratings',['All respondents, country adjusted','Heterosexual respondents','Equal country weight, all','Equal country weight, heterosexual','Country, age and sociosexuality'],[-.246641,-.204831,-.173640,-.170259,-.009791],[-.375433,-.344604,-.335079,-.347758,-.140967],[-.117850,-.065058,-.012201,.007240,.121386],[D,D,N,N,F],'Male minus female decline (rating points; 95% CI)',(-.45,.20))
forest('04_marriage_comparison.png','Marriage status and lifetime partner counts',['Date Psychology: women','Date Psychology: men','Date Psychology: men − women','NSFG: women','NSFG: men','NSFG: men − women'],[3.33,5.58,2.25,3.333229,4.080408,.747179],[-3.97,1.08,-6.29,2.702780,2.931692,-.431973],[10.63,10.08,10.79,3.963679,5.229124,1.926330],[F,M,N,F,M,N],'Unmarried minus married mean count, or difference of gaps',(-8,12))
bars('05_nsfg_status.png','Marriage associations depend on the outcome',['Currently married','Ever married','Married or\ncohabiting'],[('Men',[-16.23,-7.75,-8.33],M),('Women',[-24.18,-14.38,-15.79],F)],'Adjusted difference (percentage points)',(-30,3))
p=pd.read_csv(R/'phase11'/'expanded_design_estimates.csv');labels=[];b=[];lo=[];hi=[];cs=[]
for e,l in [('maximum','Maximum index'),('cumulative','Cumulative index'),('rawmaxgroup','Direct reports only'),('regularmaxgroup','Regular intervals or lifetime')]:
 for s,c in [('Men',M),('Women',F)]:
  q=p[(p.exposure==e)&(p.sex==s)].iloc[0];labels.append(l+' · '+s.lower());b.append(q.rd*100);lo.append(q.lo*100);hi.append(q.hi*100);cs.append(c)
forest('06_prospective_robustness.png','Lower first-marriage attainment by age 40',labels,b,lo,hi,cs,'10+ versus 1–4 at age 20 (percentage points; 95% CI)',(-26,3))
j=pd.read_csv(R/'phase11'/'joint_ipcw_estimates.csv');q=j[(j.exposure=='maximum')&(j.horizon==240)&(~j.older)&(j.cap==1)].set_index('sex');bars('07_prospective_risks.png','First marriage by age 40: standardized estimates',['Men','Women'],[('Maximum index 1–4',[100*q.loc[s,'low'] for s in ['Men','Women']],M),('Maximum index 10+',[100*q.loc[s,'high'] for s in ['Men','Women']],F)],'Estimated probability (%)',(0,85),True)
a=pd.read_csv(R/'paired_attenuation.csv');a=a[a.sex.isin(['Men','Women'])].copy();a['label']=a.exposure.map({'maximum_v2_group':'Maximum','cumulative_v2_group':'Cumulative'})+' · '+a.sex.str.lower();forest('08_paired_attenuation.png','Updating the index changes its predictive association',a.label.tolist(),(a.attenuation*100).tolist(),(a.ci95_lower*100).tolist(),(a.ci95_upper*100).tolist(),[M if s=='Men' else F for s in a.sex],'Attenuation from updating at age 25 (percentage points)',(-9,24))
# timing comparison from executed coefficients
x=pd.read_csv(R/'phase9'/'paired_timing.csv');x=x[(x.model=='expanded')&(x.term.isin(['0','2','3','4','5+']))];fig,axes=plt.subplots(2,1,figsize=(8.3,6.4),sharex=True);cats=['0','2','3','4','5+'];xx=np.arange(5)
for ax,s in zip(axes,['Men','Women']):
 t=x[x.sex==s].set_index('term').loc[cats];ax.bar(xx-.18,np.exp(t.current),.34,color=N,label='Contemporaneous');ax.bar(xx+.18,np.exp(t.lag),.34,color=M,label='Strictly lagged');ax.axhline(1,color='#8392a0',ls='--',lw=1);ax.grid(axis='y',alpha=.16);ax.set_ylim(0,1.08);ax.set_title(s,loc='left',fontsize=13);ax.set_ylabel('Hazard ratio');ax.set_xticks(xx,cats);ax.tick_params(axis='x',length=0)
axes[0].legend(loc='upper center',bbox_to_anchor=(.5,1.24),ncol=2,frameon=False);axes[-1].set_xlabel('Recent partners (reference: exactly one)');fig.suptitle('Prospective measurement attenuates the recent-partner association',x=.1,ha='left',fontsize=14,fontweight='bold');fig.subplots_adjust(top=.85,hspace=.35);save(fig,'09_timing.png')
