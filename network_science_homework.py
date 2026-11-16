"""Run the archived network science experiments and optional figures."""

import argparse
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
SECTIONS = (
    "all",
    "download",
    "er",
    "power",
    "ws",
    "real",
    "community",
    "tail",
    "validate",
)


def run(script, *args):
    subprocess.run(
        [sys.executable, str(BASE / "scripts" / script), *args], cwd=BASE, check=True
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--section", choices=SECTIONS, default="all")
    parser.add_argument("--figures-only", action="store_true", help="Use existing CSVs")
    args = parser.parse_args()
    if args.figures_only:
        run("run_experiments.py", "--section", "download")
        run("make_figures.py")
    else:
        run("run_experiments.py", "--section", args.section)
        if args.section == "all":
            run("make_figures.py")


if __name__ == "__main__":
    main()
