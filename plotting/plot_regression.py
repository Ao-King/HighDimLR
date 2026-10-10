"""Reproduce the Case I and Case II threshold-error figures from saved results."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("Landweber", "Showalter", "HBF", "ARk", "SOAR", "Nesterov", "FAR", "Ridge", "SC", "Lasso")
COLORS = dict(Landweber="#0072B2", Showalter="#E69F00", HBF="#009E73", ARk="#8A5AAB",
              SOAR="#9D6B15", Nesterov="#00A7AD", FAR="#737373", Ridge="#222222", SC="#8C8F26", Lasso="#CC4778")
STYLES = dict(Landweber="-", Showalter=(0, (4, 1.5)), HBF="-", ARk=(0, (2.5, 1.3)),
              SOAR="-", Nesterov=(0, (5, 1.2, 1, 1.2)), FAR=(0, (5, 2)), Ridge="-", SC=(0, (1, 1.2)), Lasso="-")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_case(directory, case):
    summary = read_json(directory / "summary.json")
    if summary["status"] != "complete":
        raise ValueError(f"Incomplete results: {directory}")
    audit_path = directory / "audit.json"
    if not audit_path.exists():
        audit_path = directory / "saved_results_audit.json"
    if read_json(audit_path)["status"] != "passed":
        raise ValueError(f"Results did not pass their audit: {directory}")
    entries = {}
    grid = .001 * np.arange(1, 1001)
    maximum = 0.
    for row in summary["rows"]:
        method = row["method"]
        if method not in METHODS:
            continue
        noise = row["distribution"]
        with np.load(directory / f"{noise.lower()}_{method.lower()}.npz", allow_pickle=False) as saved:
            scan = saved["threshold_scan"].copy()
            raw, fitted, thresholds, beta = (saved[key].copy() for key in
                                            ("raw_estimates", "estimates", "thresholds", "beta"))
        branches = 1 if method == "Lasso" else 2
        if scan.shape != (1000, 3) or not np.isfinite(scan).all() or not np.array_equal(scan[:, 0], grid):
            raise ValueError(f"Invalid threshold grid for {case}/{noise}/{method}")
        if raw.shape != fitted.shape or raw.shape != (summary["p"], branches):
            raise ValueError(f"Invalid estimates for {case}/{noise}/{method}")
        for j, branch in enumerate(("I", "II")[:branches]):
            curve = scan[:, j + 1]
            indices = np.flatnonzero(curve == curve.min())
            selected = (np.percentile(indices, row["percentile"])
                        if len(indices) > 1 and np.all(np.diff(indices) == 1) else indices[0])
            expected = float(np.interp(selected, np.arange(1000), grid))
            thresholded = np.where(np.abs(raw[:, j]) <= expected, 0., raw[:, j])
            if not (np.isclose(row[f"threshold_{branch}"], expected, rtol=0, atol=1e-14)
                    and np.isclose(thresholds[j], expected, rtol=0, atol=1e-14)
                    and np.isclose(row[f"l2_{branch}"], curve.min(), rtol=1e-12, atol=1e-12)
                    and np.array_equal(thresholded, fitted[:, j])
                    and np.isclose(np.linalg.norm(thresholded - beta), row[f"l2_{branch}"], rtol=1e-12, atol=1e-12)):
                raise ValueError(f"Threshold marker and saved estimates disagree for {case}/{noise}/{method}/{branch}")
        key = noise, method
        if key in entries:
            raise ValueError(f"Duplicate saved result: {key}")
        entries[key] = row, scan
        maximum = max(maximum, float(scan[:, 1:branches + 1].max()))
    expected_keys = {(noise, method) for noise in ("Normal", "Laplace") for method in METHODS}
    if set(entries) != expected_keys:
        raise ValueError(f"The complete set of curves is required for Case {case}")
    return summary, entries, (0., float(np.ceil(maximum * 1.06 / .2) * .2))


def draw_case(summary, entries, limits, case, output_directory):
    width, height = 7.2, 4.3
    fig = plt.figure(figsize=(width, height))
    left, plots_right, col_gap = .60, 5.91, .39
    bottom, top, row_gap = .49, height - .55, .40
    panel_width = (plots_right - left - col_gap) / 2
    panel_height = (top - bottom - row_gap) / 2
    axes, labels = [], []
    for row_index, noise in enumerate(("Normal", "Laplace")):
        y = top - panel_height if row_index == 0 else bottom
        for column, branch in enumerate(("I", "II")):
            letter = chr(ord("a") + 2 * row_index + column)
            x = left + column * (panel_width + col_gap)
            ax = fig.add_axes([x / width, y / height, panel_width / width, panel_height / height])
            description = "ordinary" if column == 0 else "debiased"
            labels.append(ax.set_title(f"({letter}) {noise}: {description}", fontsize=9.7, pad=5))
            for method in METHODS:
                if column == 1 and method == "Lasso":
                    continue
                row, scan = entries[noise, method]
                ax.plot(scan[:, 0], scan[:, column + 1], color=COLORS[method], linestyle=STYLES[method],
                        linewidth=1.15, zorder=2, solid_capstyle="round")
                ax.plot(row[f"threshold_{branch}"], row[f"l2_{branch}"], marker="o", markersize=4.2,
                        markeredgewidth=.9, markerfacecolor="white", markeredgecolor=COLORS[method],
                        linestyle="none", zorder=4)
                if not (np.array_equal(ax.lines[-2].get_xdata(), scan[:, 0])
                        and np.array_equal(ax.lines[-2].get_ydata(), scan[:, column + 1])):
                    raise ValueError("Plotted data differ from the saved curve")
            ax.set(xlim=(0, 1), ylim=limits)
            ax.xaxis.set_major_locator(MultipleLocator(.2))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.grid(axis="y", color="#BDBDBD", alpha=.35, linewidth=.45, zorder=0)
            ax.tick_params(direction="out", length=3, width=.6, pad=3, labelsize=8.7, labelleft=(column == 0))
            ax.spines[["top", "right"]].set_visible(False)
            formula = r"$\|\hat{\theta}-\beta\|_2$" if column == 0 else r"$\|\tilde{\theta}-\beta\|_2$"
            labels.append(ax.set_ylabel(formula, fontsize=9.7, labelpad=6))
            if row_index == 1:
                labels.append(ax.set_xlabel(r"Threshold $b_n$", fontsize=9.7, labelpad=4))
            else:
                ax.tick_params(labelbottom=False)
            axes.append(ax)
    dimensions = (rf"$n=p={summary['n']}$" if summary["n"] == summary["p"] else
                  rf"$n={summary['n']},\ p={summary['p']}$")
    labels.append(fig.text((left + plots_right) / (2 * width), .981,
                          f"Case {case}: threshold-error curves ({dimensions})",
                          ha="center", va="top", fontsize=11))
    handles = [Line2D([], [], color=COLORS[m], linestyle=STYLES[m], linewidth=1.5,
                      label=r"AR$^\kappa$" if m == "ARk" else m) for m in METHODS]
    legend = fig.legend(handles=handles, loc="center left",
                        bbox_to_anchor=(6.075 / width, (top + bottom) / (2 * height)),
                        ncol=1, frameon=False, fontsize=9.2, handlelength=1.95, handletextpad=.55,
                        labelspacing=.65, borderaxespad=0, title="Method", title_fontsize=9.5)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    if not (box.x0 > axes[1].get_position().x1 and box.x1 < 1 and box.y0 > 0 and box.y1 < 1):
        raise ValueError("The legend extends beyond its reserved area")
    for artist in labels:
        bounds = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        if not (bounds.x0 >= 0 and bounds.x1 <= 1 and bounds.y0 >= 0 and bounds.y1 <= 1):
            raise ValueError(f"A label extends beyond the figure: {artist.get_text()}")
    stem = "Fig1_CaseI" if case == "I" else "FigD1_CaseII"
    fig.savefig(output_directory / f"{stem}.pdf", metadata={"Title": f"Case {case}: threshold-error curves"})
    fig.savefig(output_directory / f"{stem}.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/regression")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/figures")
    parser.add_argument("--case", choices=("I", "II", "all"), default="all")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 10,
                         "axes.linewidth": .6, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.facecolor": "white", "path.simplify": False})
    for case in (("I", "II") if args.case == "all" else (args.case,)):
        summary, entries, limits = load_case(args.data_dir / f"case{case}" / "point", case)
        draw_case(summary, entries, limits, case, args.output_dir)
    print(args.output_dir.resolve())


if __name__ == "__main__":
    main()
