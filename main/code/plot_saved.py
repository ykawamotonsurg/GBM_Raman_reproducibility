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

def fucoidan_annotations(summary, tests):
    """Format saved effects and primary-test P values; do not recompute them."""
    effects = summary.set_index('Readout', verify_integrity=True)
    primary = tests.loc[tests['Test'] == 'Exact signed-rank (primary)'].set_index('Band (cm⁻¹)', verify_integrity=True)
    titles, pvalues = {}, {}
    for band in (1156, 1514):
        column = f'{band} cm⁻¹'
        lower = effects.loc['Pairs with a lower signal', column]
        reduction = effects.loc['Reduction from the ratio of condition means²', column]
        if pd.isna(lower) or pd.isna(reduction):
            raise ValueError(f'Missing saved effect annotation for {band}')
        titles[band] = f'{lower} pairs lower; mean {reduction} lower'
        p = float(primary.loc[band, 'BH-adjusted P'])
        if not np.isfinite(p) or not 0 <= p <= 1:
            raise ValueError(f'Invalid saved primary BH-adjusted P for {band}')
        pvalues[band] = p
    if pvalues[1156] == pvalues[1514]:
        heading = f'BH-adjusted P = {pvalues[1156]:.5g} for both bands'
    else:
        heading = '; '.join(f'{band} cm⁻¹: BH-adjusted P = {pvalues[band]:.5g}' for band in (1156, 1514))
    return titles, heading

def run(root,output):
    root=Path(root);out=Path(output)/'numeric';out.mkdir(parents=True,exist_ok=True)
    def read(name,source=False):return pd.read_csv(resolve(root,('data/source_data/' if source else 'data/processed/')+name),float_precision='round_trip')
    records=[]
    def save(fig,name,inputs,note='Fixed saved numerical values'):
        fig.savefig(out/(name+'.png'),dpi=240,bbox_inches='tight');plt.close(fig)
        records.append({'panel':name,'level':'A','inputs':inputs,'output':'numeric/'+name+'.png','operation':note})
    def axnew(w=6,h=4):return plt.subplots(figsize=(w,h),layout='constrained')
    def regx(ax):ax.set_xticks(range(4),REG);ax.set_xlabel('Region')
    bags=read('tissue_bags.csv')
    fig,ax=axnew()
    vals=[np.log(np.clip(bags[bags.region==r].bag_score.to_numpy(),1e-6,1-1e-6)/(1-np.clip(bags[bags.region==r].bag_score.to_numpy(),1e-6,1-1e-6))) for r in REG]
    ax.boxplot(vals,positions=range(4),showfliers=False,widths=.5)
    for j,v in enumerate(vals):ax.scatter(np.full(len(v),j)+np.linspace(-.08,.08,len(v)),v,c='#777777',s=14,alpha=.8)
    regx(ax);ax.set_ylabel('Tissue Raman score (logit)')
    save(fig,'Fig_1A',['main/data/tissue_bags.csv'],'Fixed probabilities clipped to [1e-6,1-1e-6], logit display only; box quartiles and deterministic display offsets')
    roc=read('tissue_roc.csv');stat=read('tissue_roc_statistics.csv')
    for comp,name in [('IM_ED','Fig_1B')]:
        z=roc[roc.comparison==comp];s=stat[(stat.comparison==comp)&(stat.bootstrap_unit=='class_stratified_bag')].iloc[0]
        fig,ax=axnew(4.8,4.8);ax.plot(z.FPR,z.TPR,c='black',lw=1.5);ax.plot([0,1],[0,1],'--',c='gray')
        ax.set(xlim=(0,1),ylim=(0,1),xlabel='1 - Specificity',ylabel='Sensitivity')
        ax.text(.97,.07,f'AUC = {s.AUC:.3f}\n95% CI {s.CI_low:.3f}–{s.CI_high:.3f}\nClass-stratified bag bootstrap',transform=ax.transAxes,ha='right',fontsize=9)
        save(fig,name,['shared/data/tissue_roc.csv','shared/data/tissue_roc_statistics.csv'],'Saved ROC coordinates and saved AUC/intervals; no bootstrap or classifier run')
    z=read('tissue_spectra.csv');fig,ax=axnew(7,4)
    for r in REG:
        q=z[z.area4==r];ax.plot(q.raman_shift,q.mean_intensity,label=r,color=PEAKCOL[r],lw=1)
    for w in [1002,1156,1440,1514,1654]:ax.axvline(w,color='gray',ls=':',lw=.65);ax.text(w,.225,str(w),rotation=90,ha='center',fontsize=8)
    ax.set(xlim=(900,1800),ylim=(0,.24),xlabel=r'Raman shift (cm$^{-1}$)',ylabel='Relative intensity');ax.legend(title='Region',bbox_to_anchor=(1,1))
    save(fig,'Fig_1C',['main/data/tissue_spectra.csv'])
    peaks=read('tissue_peaks.csv');ann=read('tissue_peak_annotations.csv')
    for pp,name in [([1156,1514],'Fig_1D')]:
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
    sp=read('cell_mean_spectra.csv');fig,ax=axnew(8,3.8)
    for group,col in [('Carotenoid-high','#00B9C4'),('Carotenoid-low','#F8766D')]:
        z=sp[sp.group==group];ax.plot(z.wavenumber,z['mean'],c=col,label=group)
    for w in [1156,1514]:ax.axvline(w,color='gray',ls='--',lw=.8);ax.text(w,140,str(w),ha='center',fontsize=9)
    ax.set(xlabel=r'Raman shift (cm$^{-1}$)',ylabel='Intensity',xlim=(900,1800),ylim=(0,150));ax.legend(bbox_to_anchor=(1,1),frameon=False)
    save(fig,'Fig_2A',['main/data/cell_mean_spectra.csv'])
    cl=read('cells.csv')
    for mode,name in [('region','Fig_2B'),('high','Fig_2C')]:
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
    z=read('Singlecell_CaseRegion.csv',True);fig,ax=axnew(4.8,4.6)
    for _,g in z.groupby('patient_id'):
        g=g.assign(x=g.region.map({r:i for i,r in enumerate(REG)})).sort_values('x')
        if len(g)>1:ax.plot(g.x,g.frac_caro_high*100,color='#dddddd',lw=.8)
    for j,r in enumerate(REG):
        q=z[z.region==r];ax.scatter(np.full(len(q),j),q.frac_caro_high*100,c=COL[r],s=24)
    regx(ax);ax.set_ylim(-2,82);ax.set_ylabel('Carotenoid-high fraction (%)')
    save(fig,'Fig_2D',['main/data/Singlecell_CaseRegion.csv'],'38 saved summaries; lines retain patient pairing')
    # Fig.3B squares read directly from the approved counts table.
    z=read('Table_S6.csv',True);im=Image.new('RGB',(1210,1530),'white');dr=ImageDraw.Draw(im)
    normal=findfont(FontProperties(family='DejaVu Sans'));bold=findfont(FontProperties(family='DejaVu Sans',weight='bold'))
    def text(t,x,y,size=40,b=False):dr.text((x,y),str(t),fill='black',font=ImageFont.truetype(bold if b else normal,size))
    text('B',10,5,66,True);text('Same-cell correspondence',95,35);text('MSR1+',235,115,56);text('MSR1−',755,115,56)
    for j,r in z.iterrows():
        y=210+j*395;text(r['Case'],25,y+88,48)
        for k,(hc,lc) in enumerate([('High MSR1+','Low MSR1+'),('High MSR1−','Low MSR1−')]):
            h=int(r[hc]);n=h+int(r[lc]);x=205+k*510;text(f'{h}/{n} high' if n else 'No evaluable cells',x,y,50 if n else 43)
            for v in range(n):
                xx=x+v%6*66;yy=y+78+v//6*66;dr.rectangle((xx,yy,xx+55,yy+55),fill='black' if v<h else 'white',outline='black',width=3)
    text('Each square = one evaluable cell',190,1410);dr.rectangle((190,1470,225,1505),fill='black');text('high',242,1468,43);dr.rectangle((525,1470,560,1505),outline='black',width=3);text('low',575,1468,43)
    im.save(out/'Fig_3B.png');records.append({'panel':'Fig_3B','level':'A','inputs':['main/data/Table_S6.csv'],'output':'numeric/Fig_3B.png','operation':'55 squares from exact stored counts; no test'})
    z=read('Culture_AllConditions.csv',True);fig=plt.figure(figsize=(10.8,3.6));axs=[fig.add_axes([.09,.24,.39,.60]),fig.add_axes([.56,.24,.39,.60])]
    for ax,model,title in zip(axs,['THP1','U87'],['THP-1-derived M2-like','U-87 MG']):
        for k,(cond,col) in enumerate([('Control','#bdbdbd'),('oxLDLcar','#f4a259')]):
            q=z[(z.model==model)&(z.Condition==cond)].sort_values('peak');x=np.array([0,1])+(k-.5)*.25
            ax.bar(x,q['mean'],width=.24,yerr=q.sd,color=col,edgecolor='#666666',lw=.8,capsize=4,label='Control' if k==0 else r'oxLDL$^{car}$')
            for xi,(_,r) in zip(x,q.iterrows()):ax.text(xi,r['mean']+r.sd+1.2,'n='+str(r['n']),ha='center',fontsize=10)
        ax.set_xticks([0,1],['1156','1514']);ax.set_ylim(0,75);ax.set_xlim(-.5,1.5);ax.set_yticks([0,20,40,60]);ax.set_title(title,fontsize=14,fontweight='bold',pad=12);ax.set_xlabel(r'Raman shift (cm$^{-1}$)')
    axs[0].set_ylabel('Peak intensity');fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.01))
    # Match the approved display build's fixed pixel canvas for the frozen assembly.
    fig.savefig(out/'Fig_3C_assembly.png',dpi=350,facecolor='white')
    save(fig,'Fig_3C',['main/data/Culture_AllConditions.csv'],'Saved mean, SD and n; control/loading selected, all conditions retained in source data')
    z=read('Fucoidan_Pairs.csv',True);fig,axs=plt.subplots(1,2,figsize=(10,4.7),layout='constrained');dc={'Day1':'#1f77b4','Day2':'#ff7f0e','Day3':'#2ca02c'}
    titles, heading = fucoidan_annotations(read('Table_S7A.csv',True), read('Table_S7B.csv',True))
    for ax,pk in zip(axs,[1156,1514]):
        title=titles[pk]
        q=z[z.peak==pk]
        for _,r in q.iterrows():ax.plot([0,1],[r.mean_intensity_caro,r.mean_intensity_fuco],'-o',color=dc[r.experiment_id],ms=4,lw=1)
        ax.set_xticks([0,1],[r'oxLDL$^{car}$',r'oxLDL$^{car}$'+'\n+ fucoidan']);ax.set(xlim=(-.2,1.2),ylim=(0,75),ylabel='Peak intensity');ax.set_title(f'{pk} cm$^{{-1}}$\n{title}',fontsize=11)
    handles=[plt.Line2D([0],[0],color=v,marker='o',label=k) for k,v in dc.items()]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False);fig.suptitle(heading,fontsize=10)
    save(fig,'Fig_3D',['main/data/Fucoidan_Pairs.csv','main/data/Table_S7A.csv','main/data/Table_S7B.csv'],'6 matched pairs across 3 saved day blocks; existing effect/P annotations only')
    xp=read('spatial_program_scores.csv');progs=list(xp.program.unique());labels=dict(zip(xp.program,xp.program_label))
    # Preserve the source program order used in the approved figure.
    assert set(progs)==set(PROGRAMS)
    progs=PROGRAMS
    def program_ax(ax,z,niches):
        for case,g in z.groupby('case',sort=True):
            g=g.set_index('niche').loc[niches].reset_index();xx=np.arange(len(niches));ax.plot(xx,g.score,c='#555555',lw=.9);ax.scatter(xx,g.score,s=10+g.n_cells.to_numpy()/220,c='#333333')
        ax.set_xticks(range(len(niches)),[{'I':'Interstitial','R':'Rim','P':'Perivasc.'}[v] for v in niches],fontsize=8);ax.grid(color='#eeeeee');ax.set_axisbelow(True)
    fig,axs=plt.subplots(2,4,figsize=(13,8),layout='constrained')
    for ax,pr in zip(axs.flat,progs):program_ax(ax,xp[(xp.program==pr)&(xp.group=='MSR1+')],['I','R']);ax.set_title(labels[pr],fontsize=10);ax.set_ylabel('Program score',fontsize=9)
    axs.flat[-1].axis('off');handles=[axs.flat[-1].scatter([],[],s=10+n/220,c='#333333',label=str(n)) for n in [2000,3000,4000,5000]];axs.flat[-1].legend(handles=handles,title='Cells (n)',loc='center',frameon=False)
    save(fig,'Fig_4B',['shared/data/spatial_program_scores.csv'],'All 7 programs; 3 fixed case-level lines; point size uses saved cell counts')
    return records
