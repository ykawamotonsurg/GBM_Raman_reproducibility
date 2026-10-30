"""Summarize recorded cell labels and culture units; no classification or testing."""
from pathlib import Path
from collections import defaultdict
import argparse, csv, json, math, statistics

ROOT=Path(__file__).resolve().parents[2]

def read(rel):
    with (ROOT/rel).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def write(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out.resolve()
    if out.exists():p.error('Choose a new output directory.')
    if any(out==ROOT/d or ROOT/d in out.parents for d in ('main','shared','supplement','metadata')):
        p.error('Output must be outside the supplied data and code folders.')
    labels=read('main/data/samecell_labels.csv');groups=defaultdict(list)
    assert len({r['samecell_id'] for r in labels})==len(labels)
    for r in labels:
        assert r['carotenoid_high'] in ('True','False')
        assert r['MSR1_pos_p99'] in ('','0.0','1.0')
        groups[r['display_case']].append(r)
    counts=[]
    for case,rs in sorted(groups.items()):
        row={'Case':case}
        for high,state,col in [(True,'1.0','High MSR1+'),(True,'0.0','High MSR1−'),
                               (False,'1.0','Low MSR1+'),(False,'0.0','Low MSR1−')]:
            row[col]=sum(r['carotenoid_high']==str(high) and r['MSR1_pos_p99']==state for r in rs)
        row['Evaluable cells']=sum(r['MSR1_pos_p99']!='' for r in rs);counts.append(row)
    expected=read('main/data/Table_S6.csv')
    assert [{k:str(v) for k,v in r.items()} for r in counts]==expected
    units=read('main/data/culture_unit_means.csv');culture=defaultdict(list)
    assert len({(r['unit_id'],r['Condition'],r['peak']) for r in units})==len(units)
    for r in units: culture[(r['model'],r['Condition'],r['peak'])].append(float(r['mean_intensity']))
    summaries=[]
    for r in read('main/data/Culture_AllConditions.csv'):
        vs=culture[(r['model'],r['Condition'],r['peak'])]
        mean=statistics.mean(vs);sd=statistics.stdev(vs)
        assert len(vs)==int(r['n'])
        assert math.isclose(mean,float(r['mean']),rel_tol=1e-12,abs_tol=1e-12)
        assert math.isclose(sd,float(r['sd']),rel_tol=1e-12,abs_tol=1e-12)
        summaries.append(dict(model=r['model'],Condition=r['Condition'],peak=r['peak'],mean=mean,sd=sd,n=len(vs)))
    cells=read('shared/data/cells.csv');cellgroups=defaultdict(list)
    for r in cells:cellgroups[(r['patient_id'],r['region'])].append(r)
    fractions=[]
    for r in read('main/data/Singlecell_CaseRegion.csv'):
        rs=cellgroups[(r['patient_id'],r['region'])];k=sum(x['carotenoid_high']=='True' for x in rs)
        assert len(rs)==int(r['n_cells']) and k==int(r['n_caro_high'])
        assert math.isclose(k/len(rs),float(r['frac_caro_high']),rel_tol=1e-12,abs_tol=1e-12)
        fractions.append(dict(patient_id=r['patient_id'],region=r['region'],n_cells=len(rs),n_caro_high=k,frac_caro_high=k/len(rs)))
    out.mkdir(parents=True)
    write(out/'Table_S6.csv',counts);write(out/'Culture_AllConditions.csv',summaries)
    write(out/'Singlecell_CaseRegion_counts.csv',fractions)
    report=dict(samecell_recorded=len(labels),samecell_evaluable=sum(r['Evaluable cells'] for r in counts),
                culture_groups=len(summaries),case_region_groups=len(fractions),
                stored_counts_and_summaries_match=True,new_classification=False,new_statistical_tests=False)
    (out/'summary_check.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report))

if __name__=='__main__':main()
