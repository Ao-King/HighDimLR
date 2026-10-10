"""Export the published or recomputed Table 3 summary to CSV and LaTeX."""
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
    parser.add_argument('--summary', type=Path, default=ROOT / 'data/regression/table3_summary.json')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/regression/table3-export')
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding='utf-8'))
    if 'methods' in summary:
        order = ('SC', 'Ridge', 'Landweber', 'Showalter', 'ARk', 'SOAR', 'Nesterov', 'HBF', 'FAR')
        rows = [dict(method=method, percentile=PERCENTILES[method],
                     **summary['methods'][method]['recomputed']) for method in order]
    else:
        rows = summary['table3']['rows'] if 'table3' in summary else summary['rows']
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / 'table3.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    tex = [r'\begin{tabular}{l' + 'r' * len(rows) + '}', r'\hline',
           ' & ' + ' & '.join(r['method'] for r in rows) + r' \\']
    labels = dict(percentile=r'$b_n^\circ$', coverage_I='Coverage I', coverage_II='Coverage II',
        l2_I=r'Mean $E_2$ I', l2_II=r'Mean $E_2$ II', linf_I=r'Mean $E_\infty$ I',
        linf_II=r'Mean $E_\infty$ II', variance_error_I='Variance error I', variance_error_II='Variance error II')
    for field in FIELDS[1:]:
        values = [f'{100*r[field]:.1f}\\%' if field.startswith('coverage') else
                  f'{r[field]:g}' if field == 'percentile' else f'{r[field]:.4f}' for r in rows]
        tex.append(labels[field] + ' & ' + ' & '.join(values) + r' \\')
    tex.extend((r'\hline', r'\end{tabular}', ''))
    (args.output / 'table3.tex').write_text('\n'.join(tex), encoding='utf-8')
    print(args.output)


if __name__ == '__main__':
    main()
