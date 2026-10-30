"""Plot fixed processed values. No inference, classification or embedding fitting.

The two MRI straight lines replay the existing descriptive geom_smooth(method='lm')
display on the same saved points. They do not produce tests or confidence intervals.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"shared/code"))
from io_paths import resolve

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont

REG=['CWM','ED','IM','CE']
COL=dict(zip(REG,['#BDBDBD','#67A9CF','#F4A259','#EF4135']))
PEAKCOL=dict(zip(REG,['#F8766D','#7CAE00','#00BFC4','#C77CFF']))
PROGRAMS=['HypoxiaGlyco','Stress','Opsonin_sensing','Lysosomal_degradation','Uptake_receptor','Lipid_handling','Efflux_export']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})

def run(root,output):
    root=Path(root);out=Path(output)/'numeric';out.mkdir(parents=True,exist_ok=True)
    def read(name,source=False):return pd.read_csv(resolve(root,('data/source_data/' if source else 'data/processed/')+name),float_precision='round_trip')
    records=[]
    def save(fig,name,inputs,note='Fixed saved numerical values'):
        fig.savefig(out/(name+'.png'),dpi=240,bbox_inches='tight');plt.close(fig)
        records.append({'panel':name,'level':'A','inputs':inputs,'output':'numeric/'+name+'.png','operation':note})
    def axnew(w=6,h=4):return plt.subplots(figsize=(w,h),layout='constrained')
    def regx(ax):ax.set_xticks(range(4),REG);ax.set_xlabel('Region')
    roc=read('tissue_roc.csv');stat=read('tissue_roc_statistics.csv')
    for comp,name in [('CE_CWM','Fig_S2A')]:
        z=roc[roc.comparison==comp];s=stat[(stat.comparison==comp)&(stat.bootstrap_unit=='class_stratified_bag')].iloc[0]
        fig,ax=axnew(4.8,4.8);ax.plot(z.FPR,z.TPR,c='black',lw=1.5);ax.plot([0,1],[0,1],'--',c='gray')
        ax.set(xlim=(0,1),ylim=(0,1),xlabel='1 - Specificity',ylabel='Sensitivity')
        ax.text(.97,.07,f'AUC = {s.AUC:.3f}\n95% CI {s.CI_low:.3f}–{s.CI_high:.3f}\nClass-stratified bag bootstrap',transform=ax.transAxes,ha='right',fontsize=9)
        save(fig,name,['shared/data/tissue_roc.csv','shared/data/tissue_roc_statistics.csv'],'Saved ROC coordinates and saved AUC/intervals; no bootstrap or classifier run')
    peaks=read('tissue_peaks.csv');ann=read('tissue_peak_annotations.csv')
    for pp,name in [([1440,1654],'Fig_S2C')]:
        fig,axs=plt.subplots(1,2,figsize=(8,4),layout='constrained')
        for ax,pk in zip(axs,pp):
            z=peaks[peaks.peak_cm1==pk]
            for _,g in z.groupby('patient_id'):
                g=g.assign(x=g.region.map({r:i for i,r in enumerate(REG)})).sort_values('x')
                if len(g)>1:ax.plot(g.x,g.value,color='#dddddd',lw=.7,zorder=1)
            ax.boxplot([z[z.region==r].value for r in REG],positions=range(4),showfliers=False,widths=.5)
            for j,r in enumerate(REG):ax.scatter(np.full(sum(z.region==r),j),z[z.region==r].value,c=PEAKCOL[r],s=12,zorder=3)
            a=ann[ann.peak_cm1==pk].iloc[0];ax.plot([1,1,2,2],[a.y-a.tip,a.y,a.y,a.y-a.tip],c='black',lw=.7);ax.text(1.5,a.ytext,a.label,fontsize=7,ha='center')
            ax.set_ylim(0,max(z.value.max()*1.2,a.ytext*1.07));ax.set_title(str(pk),fontsize=12);regx(ax)
        axs[0].set_ylabel('Relative intensity');save(fig,name,['shared/data/tissue_peaks.csv','shared/data/tissue_peak_annotations.csv'],'57 patient-region values per band, 39 patient IDs; existing Holm annotations loaded, no tests')
    cl=read('cells.csv')
    for mode,name in [('cluster','Fig_S3A')]:
        fig,ax=axnew(5.4,4.7)
        if mode=='high':
            ax.scatter(cl.UMAP1,cl.UMAP2,c='#dddddd',s=3)
            z=cl[cl.carotenoid_high];ax.scatter(z.UMAP1,z.UMAP2,c='#009C68',s=4)
        else:
            groups=REG if mode=='region' else sorted(cl.cluster.unique())
            palette=COL if mode=='region' else dict(zip(groups,['#619CFF','#F8766D','#00BA38','#C77CFF','#00BFC4','#B79F00','#FF61CC','#8B6D5C']))
            for group in groups:
                z=cl[cl[mode]==group];ax.scatter(z.UMAP1,z.UMAP2,c=palette[group],s=3,alpha=.85,label=str(group))
            ax.legend(title='Region' if mode=='region' else 'Cluster',bbox_to_anchor=(1,1),frameon=False,markerscale=2)
        ax.set_aspect('equal');ax.axis('off')
        if mode=='cluster':ax.set_title('Raman-defined cluster')
        save(fig,name,['shared/data/cells.csv'],'All 1,499 saved coordinates and original labels; no clustering or UMAP calculation')
    xp=read('spatial_program_scores.csv');progs=list(xp.program.unique());labels=dict(zip(xp.program,xp.program_label))
    # Preserve the source program order used in the approved figure.
    assert set(progs)==set(PROGRAMS)
    progs=PROGRAMS
    def program_ax(ax,z,niches):
        for case,g in z.groupby('case',sort=True):
            g=g.set_index('niche').loc[niches].reset_index();xx=np.arange(len(niches));ax.plot(xx,g.score,c='#555555',lw=.9);ax.scatter(xx,g.score,s=10+g.n_cells.to_numpy()/220,c='#333333')
        ax.set_xticks(range(len(niches)),[{'I':'Interstitial','R':'Rim','P':'Perivasc.'}[v] for v in niches],fontsize=8);ax.grid(color='#eeeeee');ax.set_axisbelow(True)
    z=read('Singlecell_Pooled.csv',True).set_index('region').loc[REG];fig,ax=axnew(5,4.5);x=np.arange(4)
    ax.errorbar(x,z.frac*100,yerr=[(z.frac-z.low)*100,(z.high-z.frac)*100],fmt='o',color='black',capsize=3)
    for i,(_,r) in enumerate(z.iterrows()):ax.text(i,r.high*100+1.5,f'{r.frac*100:.1f}%',ha='center')
    ax.set_xticks(x,[f'{r}\n(n={int(z.loc[r,"n"])})' for r in REG]);ax.set(ylabel='Carotenoid-high fraction (%)',ylim=(-1.5,30))
    save(fig,'Fig_S3B',['supplement/data/Singlecell_Pooled.csv'],'Saved fractions and Wilson intervals; no interval recalculation')
    z=read('cell_sensitivity.csv');fig,axs=plt.subplots(3,3,figsize=(10,8),layout='constrained')
    for i,score in enumerate(['sum','z1156','z1514']):
        for j,frac in enumerate([.1,.15,.2]):
            ax=axs[i,j];q=z[(z.score_method==score)&(z.top_fraction==frac)].set_index('region').loc[REG]
            ax.bar(range(4),q.fraction*100,color=[COL[r] for r in REG],width=.62)
            for k,v in enumerate(q.fraction):ax.text(k,100*v+1.2,f'{100*v:.1f}',ha='center',fontsize=8)
            ax.set(ylim=(0,38));regx(ax)
            if i==0:ax.set_title(f'Top {int(frac*100)}%')
            if j==0:ax.set_ylabel({'sum':'z(1156)+z(1514)','z1156':'z(1156)','z1514':'z(1514)'}[score]+'\nHigh fraction (%)')
    save(fig,'Fig_S3C',['supplement/data/cell_sensitivity.csv'],'Stored sensitivity summaries; no new cutoff or relabeling')
    z=read('tissue_cell_bridge.csv');fig,ax=axnew(7,4.5)
    for r in REG:
        q=z[z.region==r];ax.scatter(q.frac_caro_high*100,q.bag_logit,c=COL[r],label=r,s=25)
    ax.set(xlabel='Carotenoid-high fraction (%)',ylabel='Tissue Raman score (logit)',title='All regions: n = 37, Spearman rho = 0.470\nED/IM/CE: n = 33, Spearman rho = 0.368');ax.legend(title='Region',bbox_to_anchor=(1,1),frameon=False)
    save(fig,'Fig_S3D',['supplement/data/tissue_cell_bridge.csv','supplement/data/bridge_descriptive_correlations.csv'])
    z=read('spatial_hypoxia.csv').sort_values('case');fig,ax=axnew(5,4);ax.bar(z['case'],z.mean12_log_expression,color='#555555');ax.set(ylabel='Mean 12-gene log expression',ylim=(0,1.6))
    save(fig,'Fig_S4B',['supplement/data/spatial_hypoxia.csv'],'Domain-construction 12-gene saved scores, not screening HypoxiaScore')
    z=read('spatial_sensitivity.csv');pr=list(z.program_label.unique());settings=list(z.panel_label.unique());a=z.pivot(index='program_label',columns='panel_label',values='mean_effect').loc[pr,settings];ns=z.pivot(index='program_label',columns='panel_label',values='label_nsame').loc[pr,settings]
    fig,ax=axnew(10,6.7);mx=max(abs(a.to_numpy()).max(),1.2);cmap=LinearSegmentedColormap.from_list('saved_effect_display',['#2C6DB2','#FFFFFF','#E6862E']);im=ax.imshow(a,cmap=cmap,vmin=-mx,vmax=mx,aspect='auto')
    for i in range(len(pr)):
        for j in range(len(settings)):ax.text(j,i,ns.iloc[i,j],ha='center',va='center',fontsize=9)
    ax.set_yticks(range(len(pr)),pr);ax.set_xticks(range(len(settings)),[s.replace(' + ',' +\n').replace(' ','\n',1) for s in settings],fontsize=9);fig.colorbar(im,ax=ax,label='Mean difference (Rim − Interstitial)',shrink=.6)
    save(fig,'Fig_S4C',['supplement/data/spatial_sensitivity.csv'],'Stored effect means and exact sign-concordance counts; zero effects not recounted')
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,niche in zip(axs,['I','R']):
        for i,pr in enumerate(progs):
            q=xp[(xp.niche==niche)&(xp.program==pr)];neg=q[q.group=='MSR1-'].score.mean();pos=q[q.group=='MSR1+'].score.mean();ax.plot([neg,pos],[i,i],color='#555555');ax.scatter([neg,pos],[i,i],c=['#2C6DB2','#E6862E'],s=24)
        ax.set_yticks(range(7),[labels[p] for p in progs] if niche=='I' else []);ax.invert_yaxis();ax.set_xlabel('Program score');ax.set_title('Interstitial' if niche=='I' else 'Rim');ax.grid(color='#eeeeee');ax.set_axisbelow(True)
    fig.legend(handles=[plt.Line2D([0],[0],color=c,marker='o',ls='',label=n) for c,n in [('#2C6DB2','MSR1−'),('#E6862E','MSR1+')]],loc='outside lower center',ncol=2,frameon=False)
    save(fig,'Fig_S5A',['shared/data/spatial_program_scores.csv'],'Same equal-case-weight means of the 3 fixed case scores as the existing plotting script; no score fitting')
    z=read('Xenium_MSR1_Removal.csv',True);z=z[z['case'].isin(['GBM-A','GBM-B','GBM-C'])];fig,axs=plt.subplots(1,2,figsize=(9,5),layout='constrained');cc={'GBM-A':'#0072B2','GBM-B':'#D55E00','GBM-C':'#009E73'}
    for ax,niche in zip(axs,['I','R']):
        for _,r in z[z.niche==niche].iterrows():ax.plot([0,1],[r.delta_full5,r.delta_without_MSR1_4],'-o',c=cc[r['case']],label=r['case'])
        ax.axhline(0,c='gray',ls='--',lw=.7);ax.set_xticks([0,1],['Original:\n5 genes','Without MSR1:\n4 genes']);ax.set_ylim(-.05,1.55);ax.set_title('Interstitial' if niche=='I' else 'Rim')
    axs[0].set_ylabel('Case mean difference (MSR1-positive minus negative)');fig.legend(*axs[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
    save(fig,'Fig_S5B',['supplement/data/Xenium_MSR1_Removal.csv'],'Frozen full5 and four-gene contrasts; no cell selection or rescoring')
    fig,axs=plt.subplots(2,7,figsize=(18,7.5),sharey='row',layout='constrained')
    for i,group in enumerate(['MSR1+','MSR1-']):
        for j,pr in enumerate(progs):
            ax=axs[i,j];program_ax(ax,xp[(xp.program==pr)&(xp.group==group)],['P','I','R']);ax.tick_params(labelsize=7)
            if i==0:ax.set_title(labels[pr],fontsize=8)
            if j==0:ax.set_ylabel(group+'\nProgram score')
    size_handles=[axs[0,0].scatter([],[],s=10+n/220,c='#333333',label=f'{n:,}') for n in [2500,5000,7500,10000]]
    fig.legend(handles=size_handles,title='TAM cells (n)',loc='outside lower center',ncol=4,frameon=False,fontsize=8,title_fontsize=9)
    save(fig,'Fig_S6A',['shared/data/spatial_program_scores.csv'],'All 7 programs, 3 niches, both MSR1 states, all 3 primary cases; fixed scores')
    z=read('spatial_gene_shifts.csv');fig,axs=plt.subplots(1,7,figsize=(18,6),sharex=True,layout='constrained')
    for ax,pr in zip(axs,progs):
        q=z[z.program_label==labels[pr]].sort_values('delta_R_minus_I',ascending=False)
        ax.hlines(range(len(q)),0,q.delta_R_minus_I,color='#666666',lw=1);ax.scatter(q.delta_R_minus_I,range(len(q)),s=14,c='black');ax.axvline(0,c='gray',ls='--',lw=.7);ax.set_yticks(range(len(q)),q.gene,fontsize=8);ax.invert_yaxis();ax.set_title(labels[pr],fontsize=8);ax.tick_params(axis='x',labelsize=8)
    fig.supxlabel('Δ mean expression (Rim − Interstitial)');save(fig,'Fig_S6B',['supplement/data/spatial_gene_shifts.csv'],'Frozen gene-summary values only; upstream cross-case aggregation not independently established')
    agree=read('mri_agreement.csv');fig,axs=plt.subplots(2,2,figsize=(9,8.5),layout='constrained')
    for j,(metric,title) in enumerate([('tumor_volume','Tumor volume'),('necrotic_volume','Necrotic volume')]):
        z=read('mri_'+metric+'.csv');x=z.mean_xy-z.diff_xy/2;y=z.mean_xy+z.diff_xy/2;ax=axs[0,j];ax.scatter(x,y,c='black',s=20)
        lo=min(x.min(),y.min());hi=max(x.max(),y.max());ax.plot([lo,hi],[lo,hi],'k--',lw=1);coef=np.polyfit(x,y,1);xx=np.array([x.min(),x.max()]);ax.plot(xx,np.polyval(coef,xx),c='#4477ee',lw=1)
        s=agree[agree.variable==title].iloc[0];ax.text(.05,.9,f'Pearson r = {s.pearson_r:.3f}\nICC(2,1) = {s.ICC2_1:.3f} [{s.ICC2_1_lwr:.3f}, {s.ICC2_1_upr:.3f}]',transform=ax.transAxes,va='top',fontsize=8)
        ax.set(title=title,xlabel='Rater 1 volume (mm³)',ylabel='Rater 2 volume (mm³)');ax=axs[1,j];ax.scatter(z.mean_xy,z.diff_xy,c='black',s=20)
        for col,ls in [('bias','-'),('loa_upper','--'),('loa_lower','--')]:ax.axhline(z[col].iloc[0],c='black',ls=ls,lw=1)
        ax.text(.96,.68,f'Bias = {z.bias.iloc[0]:.2f}\nUpper LoA = {z.loa_upper.iloc[0]:.2f}\nLower LoA = {z.loa_lower.iloc[0]:.2f}',ha='right',transform=ax.transAxes,fontsize=8)
        ax.set(title='Bland–Altman: '+title.lower(),xlabel='Mean of two raters (mm³)',ylabel='Difference (rater 2 − rater 1) (mm³)')
    save(fig,'Fig_S7A',['supplement/data/mri_tumor_volume.csv','supplement/data/mri_necrotic_volume.csv','supplement/data/mri_agreement.csv'],'Invert saved means/differences algebraically; replay descriptive lm display; all ICC, Pearson, bias and LoA values loaded')
    z=read('MRI_CaseData.csv',True);fig,ax=axnew(6.2,4.6);x=z.enhance_fraction.to_numpy();y=100*z.fraction.to_numpy();ax.scatter(x,y,c='#303030',s=25);xx=np.array([x.min(),x.max()]);ax.plot(xx,np.polyval(np.polyfit(x,y,1),xx),c='#303030',lw=1)
    ax.set(xlim=(0,1.01),ylim=(-2,65),xlabel=r'Enhancing fraction, $V_{CE}/(V_{CE}+V_{NEC})$',ylabel='CE carotenoid-high fraction (%)',title='n = 18; Spearman rho = 0.649; P = 0.00354')
    save(fig,'Fig_S7B',['supplement/data/MRI_CaseData.csv','supplement/data/mri_spearman.csv'],'Replay existing descriptive lm line on fixed 18 points; loaded inference, no tests or bootstrap')
    return records
