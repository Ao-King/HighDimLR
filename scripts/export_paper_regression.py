"""Export Table 3 or its independent validation to CSV and LaTeX."""
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from paper_regression_config import PERCENTILES
FIELDS = ('method', 'percentile', 'coverage_I', 'coverage_II', 'l2_I', 'l2_II',
          'linf_I', 'linf_II', 'variance_error_I', 'variance_error_II')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', choices=('baseline', 'independent'), default='baseline')
    parser.add_argument('--summary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.summary is None:
        args.summary = ROOT / ('data/regression/table3_summary.json' if args.experiment == 'baseline'
                              else 'data/regression/independent_validation/summary.json')
    if args.output is None:
        name = 'table3-export' if args.experiment == 'baseline' else 'independent-validation-export'
        args.output = ROOT / 'outputs/regression' / name
    summary = json.loads(args.summary.read_text(encoding='utf-8'))
    if 'methods' in summary:
        order = ('SC', 'Ridge', 'Landweber', 'Showalter', 'ARk', 'SOAR', 'Nesterov', 'HBF', 'FAR')
        rows = [dict(method=method, percentile=PERCENTILES[method],
                     **summary['methods'][method]['recomputed']) for method in order]
    else:
        entries = summary['table3']['rows'] if 'table3' in summary else summary['rows']
        rows = []
        for entry in entries:
            row = dict(entry.get('table_row', entry))
            if 'coverage_mcse' in entry:
                row['coverage_mcse_I'], row['coverage_mcse_II'] = entry['coverage_mcse']
            rows.append(row)
    independent = (args.experiment == 'independent' or summary.get('purpose') == 'independent_validation'
                   or summary.get('experiment') == 'independent'
                   or summary.get('version', '').startswith('independent-'))
    fields = FIELDS + ('coverage_mcse_I', 'coverage_mcse_II') if independent else FIELDS
    if independent:
        count = int(summary['N'])
        for row in rows:
            for roman in ('I', 'II'):
                coverage = row['coverage_' + roman]
                row.setdefault('coverage_mcse_' + roman, (coverage * (1-coverage) / count)**.5)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / 'table3.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    tex = [r'\begin{tabular}{l' + 'r' * len(rows) + '}', r'\hline',
           ' & ' + ' & '.join(r['method'] for r in rows) + r' \\']
    labels = dict(percentile=r'$b_n^\circ$', coverage_I='Coverage I', coverage_II='Coverage II',
        l2_I=r'Mean $E_2$ I', l2_II=r'Mean $E_2$ II', linf_I=r'Mean $E_\infty$ I',
        linf_II=r'Mean $E_\infty$ II', variance_error_I='Variance error I', variance_error_II='Variance error II',
        coverage_mcse_I='Coverage MCSE I (pp)', coverage_mcse_II='Coverage MCSE II (pp)')
    for field in fields[1:]:
        values = [f'{100*r[field]:.3f}' if field.startswith('coverage_mcse') else
                  f'{100*r[field]:.1f}\\%' if field.startswith('coverage') else
                  f'{r[field]:g}' if field == 'percentile' else f'{r[field]:.4f}' for r in rows]
        tex.append(labels[field] + ' & ' + ' & '.join(values) + r' \\')
    tex.extend((r'\hline', r'\end{tabular}', ''))
    (args.output / 'table3.tex').write_text('\n'.join(tex), encoding='utf-8')
    print(args.output)


if __name__ == '__main__':
    main()
