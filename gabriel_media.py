"""Write Gabriel tissue media as frequency lists for Feko (Week 1 Demo 1 format).

Usage:  python gabriel_media.py [--out media] [--fmin-MHz 100] [--fmax-MHz 10000] [--n 150]

For each tissue (skin, fat, muscle) it writes, comma-delimited, no header:
    <tissue>_feko3col.csv   frequency [Hz], eps_r, sigma [S/m]  (CADFEKO 2022 "Import points": three columns)
and the two-column files:
    <tissue>_eps.csv     frequency [Hz], relative permittivity eps_r
    <tissue>_sigma.csv   frequency [Hz], conductivity sigma [S/m]
and course_frequencies_check.csv (with header) holding eps_r and sigma at 403, 915, 2450 and
5800 MHz. The four course frequencies are always among the list points.

The course media pack is the source to use in the report. These files are a cross-check or a
fallback; compare the Feko .out file values at 2450 MHz with course_frequencies_check.csv
(Week 1: within 1 %). If your Feko version wants another column order, units or a header,
adjust the writer below.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from tissue import gabriel_props

TISSUES = ["skin", "fat", "muscle"]
COURSE = [403e6, 915e6, 2450e6, 5800e6]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="media")
    ap.add_argument("--fmin-MHz", type=float, default=100.0)
    ap.add_argument("--fmax-MHz", type=float, default=10000.0)
    ap.add_argument("--n", type=int, default=150)
    args = ap.parse_args()

    f = np.geomspace(args.fmin_MHz * 1e6, args.fmax_MHz * 1e6, args.n)
    f = np.unique(np.concatenate([f, COURSE]))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for t in TISSUES:
        props = np.array([gabriel_props(t, x) for x in f])
        for col, name, fmt in ((0, "eps", "{:.8g}"), (1, "sigma", "{:.8g}")):
            with open(out / f"{t}_{name}.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                for x, v in zip(f, props[:, col]):
                    w.writerow([f"{x:.10g}", fmt.format(v)])

    # one file per tissue with all three columns, for the CADFEKO "Import points" button of the
    # frequency-list medium (columns: frequency, relative permittivity, conductivity)
    for t in TISSUES:
        props = np.array([gabriel_props(t, x) for x in f])
        with open(out / f"{t}_feko3col.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            for x, er, sg in zip(f, props[:, 0], props[:, 1]):
                w.writerow([f"{x:.10g}", f"{er:.8g}", f"{sg:.8g}"])

    with open(out / "course_frequencies_check.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tissue", "f_MHz", "eps_r", "sigma_Spm"])
        for t in TISSUES:
            for x in COURSE:
                er, sg = gabriel_props(t, x)
                w.writerow([t, f"{x/1e6:g}", f"{er:.4f}", f"{sg:.4f}"])
    print(f"Wrote {len(f)} points per file for {', '.join(TISSUES)} to {out}/")
    print(f"{'tissue':<8}{'f MHz':>8}{'eps_r':>9}{'sigma':>9}")
    for t in TISSUES:
        for x in COURSE:
            er, sg = gabriel_props(t, x)
            print(f"{t:<8}{x/1e6:>8g}{er:>9.2f}{sg:>9.3f}")


if __name__ == "__main__":
    main()
