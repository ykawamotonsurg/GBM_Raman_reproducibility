"""Draw saved final values; no training, threshold selection or statistical tests."""
from pathlib import Path
import argparse, csv, hashlib, importlib.util, json
ROOT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scope', choices=['main', 'supplement', 'all', 's2b'], default='main')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if out == ROOT or out in ROOT.parents or any(out == ROOT / d or ROOT / d in out.parents for d in ['main', 'supplement', 'shared', 'metadata']):
        parser.error('Use a new output directory outside the delivered content folders.')
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error('Output must be a new or empty directory.')
    with (ROOT / 'metadata/input_files.csv').open(encoding='utf-8-sig', newline='') as f:
        inputs = list(csv.DictReader(f))
    for row in inputs:
        if sha(ROOT / row['file']) != row['sha256']:
            raise ValueError('Input hash mismatch: ' + row['file'])
    out.mkdir(parents=True, exist_ok=True)
    groups = ['main', 'supplement'] if args.scope == 'all' else ['supplement'] if args.scope == 's2b' else [args.scope]
    report = {'scope': args.scope, 'input_hashes_checked': len(inputs),
              'scientific_cohort_audit_repeated': False, 'model_fit': False,
              'threshold_selection': False, 'upstream_analysis_reproduction': False,
              'S2B': 'saved-result drawing' if 'supplement' in groups else 'not requested', 'groups': {}}
    for group in groups:
        numeric = []
        if args.scope != 's2b':
            path = ROOT / group / 'code/plot_saved.py'
            spec = importlib.util.spec_from_file_location(group + '_plot', path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            numeric = module.run(ROOT, out / group)
        if group == 'supplement':
            spec = importlib.util.spec_from_file_location('s2b_plot', ROOT / 'supplement/code/plot_s2b.py')
            s2b = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(s2b)
            numeric.extend(s2b.run(ROOT, out / group))
        report['groups'][group] = {'numeric_panels': numeric}
    (out / 'execution_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'completed': groups, 'input_files': len(inputs), 'panels': sum(len(g['numeric_panels']) for g in report['groups'].values()), 'S2B': 'saved-result drawing' if 'supplement' in groups else 'not requested'}))

if __name__ == '__main__':
    main()
