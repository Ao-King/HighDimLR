"""Reproduce the seven PDE figures with MATLAB and saved plotting data."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def matlab_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matlab", default=shutil.which("matlab"),
                        help="MATLAB executable; defaults to matlab on PATH")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/plotting/pde")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/figures")
    args = parser.parse_args()
    if not args.matlab:
        parser.error("MATLAB was not found; specify --matlab or add it to PATH")
    expression = (
        f"addpath({matlab_string(Path(__file__).resolve().parent)});"
        f"render_pde({matlab_string(args.data_dir.resolve())},"
        f"{matlab_string(args.output_dir.resolve())});"
    )
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    subprocess.run([args.matlab, "-batch", expression], check=True, creationflags=flags)


if __name__ == "__main__":
    main()
