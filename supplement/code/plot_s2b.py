"""Redraw S2B from fixed late-period bag scores and saved classifications.

This is saved-result drawing, not model training, prediction, or threshold selection.
The main-analysis OOF scores are not inputs. No historical model is reconstructed.
"""
from pathlib import Path
from decimal import Decimal
import argparse
import csv
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter

DATA_FILES = ['supplement/data/s2b_temporal_late_bags.csv',
              'supplement/data/s2b_temporal_threshold.csv']
REGIONS = ['CWM', 'IM', 'ED', 'CE']
CLASSES = ['non_tumor', 'tumor']

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def load_saved(root):
    bags = read_csv(root / DATA_FILES[0])
    thresholds = read_csv(root / DATA_FILES[1])
    if len(thresholds) != 1:
        raise ValueError('Exactly one fixed threshold is required.')
    fixed = thresholds[0]
    if fixed['aggregation'] != 'top3mean' or fixed['positive_comparator'] != '>=':
        raise ValueError('Unexpected fixed aggregation or comparator.')
    if fixed['positive_class'] != 'tumor' or fixed['negative_class'] != 'non_tumor':
        raise ValueError('Unexpected saved class definition.')
    seen = set()
    cutoff = Decimal(fixed['threshold'])
    for row in bags:
        key = (row['patient_id'], row['region'])
        if key in seen or row['bag_id'] != '_'.join(key):
            raise ValueError('Duplicate or inconsistent patient-region bag.')
        seen.add(key)
        if row['region'] not in REGIONS or row['threshold_id'] != fixed['threshold_id']:
            raise ValueError('Unrecognized region or threshold reference.')
        score = Decimal(row['bag_score'])
        if not Decimal(0) <= score <= Decimal(1) or int(row['n_spectra']) <= 0:
            raise ValueError('Invalid saved score or spectrum count.')
        # Verify consistency only. Plot the stored class; never replace it or choose a threshold.
        expected = 'tumor' if score >= cutoff else 'non_tumor'
        if row['predicted_class'] != expected:
            raise ValueError('Saved classification disagrees with the supplied fixed threshold.')
    return bags, fixed

def run(root, output):
    root, output = Path(root), Path(output)
    bags, fixed = load_saved(root)
    out = output / 'numeric'
    out.mkdir(parents=True, exist_ok=True)
    counts = []
    fractions = np.zeros((len(REGIONS), len(CLASSES)))
    for i, region in enumerate(REGIONS):
        group = [r for r in bags if r['region'] == region]
        if not group:
            raise ValueError('A required region has no saved bags.')
        for j, label in enumerate(CLASSES):
            n = sum(r['predicted_class'] == label for r in group)
            counts.append({'region': region, 'predicted_class': label,
                           'n': n, 'total': len(group)})
            fractions[i, j] = n / len(group)
    style = {'font.family':'DejaVu Sans', 'font.size':12,
             'axes.spines.top':True, 'axes.spines.right':True,
             'savefig.facecolor':'white'}
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(8, 4.8), layout='constrained')
        cmap = LinearSegmentedColormap.from_list('s2b_saved_fraction', ['#132B43', '#56B1F7'])
        artist = ax.imshow(fractions, vmin=0, vmax=1, cmap=cmap, aspect='auto')
        for i, region in enumerate(REGIONS):
            for j, label in enumerate(CLASSES):
                row = counts[i * len(CLASSES) + j]
                frac = fractions[i, j]
                ax.text(j, i, f"{row['n']}/{row['total']}\n{frac:.0%}",
                        ha='center', va='center', fontsize=13,
                        color='white' if frac < .45 else 'black')
        ax.set_xticks(range(2), CLASSES)
        ax.set_yticks(range(4), REGIONS)
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True area')
        ax.set_xticks(np.arange(-.5, 2, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
        ax.grid(which='minor', color='white', linewidth=1.5)
        ax.tick_params(which='minor', bottom=False, left=False)
        bar = fig.colorbar(artist, ax=ax, fraction=.055, pad=.05,
                           ticks=[0, .25, .5, .75, 1])
        bar.ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        bar.set_label('Fraction')
        fig.savefig(out / 'Fig_S2B.png', dpi=240, bbox_inches='tight')
        plt.close(fig)
    with (out / 'Fig_S2B_counts.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['region', 'predicted_class', 'n', 'total'])
        writer.writeheader()
        writer.writerows(counts)
    return [{'panel':'Fig_S2B', 'level':'A', 'inputs':DATA_FILES,
             'output':'numeric/Fig_S2B.png', 'counts_output':'numeric/Fig_S2B_counts.csv',
             'operation':'Count saved late-period classes and draw fractions; fixed threshold consistency check only',
             'bags':len(bags), 'patients':len({r['patient_id'] for r in bags}),
             'aggregation':fixed['aggregation'], 'threshold':fixed['threshold'],
             'model_fit':False, 'model_predict':False, 'threshold_selection':False}]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    target = args.out.resolve()
    if target == root or target in root.parents or any(target == root / d or root / d in target.parents for d in ['main','supplement','shared','metadata']):
        parser.error('Use a new output directory outside delivered content folders.')
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        parser.error('Output must be a new or empty directory.')
    print(json.dumps(run(root, target), ensure_ascii=False))
