"""Export the reported PDE statistics as CSV and a LaTeX table."""

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("Landweber", "Showalter", "HBF", "ARk", "SOAR", "Nesterov", "FAR")


def scientific(value):
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def table_rows(methods):
    rows = [("tau", r"$\varsigma$", [methods[m]["options"]["tau"] for m in METHODS], "tau")]
    for branch, label in (("original", "I"), ("debiased", "II")):
        values = [100 * methods[m]["branches"][branch]["coverage"]["probability"] for m in METHODS]
        rows.append((f"coverage_{label}_percent", f"Coverage {label}", values, "coverage"))
    for metric, symbol in (("E2", "E_2"), ("Einf", r"E_\infty")):
        for branch, estimator in (("original", r"\hat f"), ("debiased", r"\tilde f")):
            values = [methods[m]["branches"][branch]["metrics"][metric]["mean"] for m in METHODS]
            rows.append((f"{metric}_{branch}", rf"$\overline{{{symbol}({estimator})}}$", values, "error"))
    for branch, estimator in (("original", r"\widehat{\sigma}"), ("debiased", r"\widetilde{\sigma}")):
        values = [methods[m]["branches"][branch]["metrics"]["absolute_variance_error"]["mean"] for m in METHODS]
        rows.append((f"absolute_variance_error_{branch}",
                     rf"$\overline{{\lvert {estimator}^2-\sigma^2\rvert}}$", values, "variance"))
    if not all(math.isfinite(float(value)) for row in rows for value in row[2]):
        raise ValueError("Table statistics must be finite")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "data/pde/table6_summary.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/tables")
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    methods = source["methods"]
    rows = table_rows(methods)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "table6.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["statistic", *METHODS])
        for key, _, values, _ in rows:
            writer.writerow([key, *values])
    parameters = methods[METHODS[0]]["parameters"]
    confidence = parameters["confidence"] * 100
    caption = (
        "Reconstruction and variance-estimation errors and empirical coverage in the inverse source experiment. "
        "Coverage I and II refer to $\\hat f$ and $\\tilde f$, respectively. "
        f"Overbars denote means over $N={parameters['N']}$ outer simulations, "
        f"with $B={parameters['B']}$ bootstrap replicates and nominal coverage ${confidence:g}\\%$."
    )
    lines = [r"\begin{table*}", r"\centering", rf"\caption{{{caption}}}", r"\label{CoveragePDE}",
             r"\begin{tabular}{c|c|c|c|c|c|c|c}", r"\hline",
             r" & Landweber & Showalter & HBF & $\mathrm{AR}^{\kappa}$ & SOAR & Nesterov & FAR \\ \hline"]
    for _, label, values, style in rows:
        if style == "coverage":
            cells = [f"{value:.1f}" + r"\%" for value in values]
        elif style == "variance":
            cells = [scientific(value) for value in values]
        elif style == "error":
            cells = [f"{value:.6f}" for value in values]
        else:
            cells = [f"{value:g}" for value in values]
        lines.append(label + " & " + " & ".join(cells) + r" \\")
    lines.extend([r"\hline", r"\end{tabular}", r"\end{table*}", ""])
    (args.output_dir / "table6.tex").write_text("\n".join(lines), encoding="utf-8")
    print(args.output_dir.resolve())


if __name__ == "__main__":
    main()
