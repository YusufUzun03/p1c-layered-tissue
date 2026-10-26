"""Compare the course media files (Gabriel lists) with this repo's Cole-Cole model.

Usage:
  python check_course_media.py --muscle Muscle.csv --skin Skin_dry.csv --fat Fat_not_infiltrated.csv
  python check_course_media.py --muscle Muscle.csv --write-3col out_dir     # also write Feko 3-column files

Files are CSV with columns: frequency [Hz], relative permittivity, conductivity [S/m] (a header line
and a fourth loss-tangent column are ignored). The course files are not part of this repo.

Prints, per tissue, the largest relative deviation of eps_r and sigma between the file and the model
at the file's own frequencies, then what Feko's linear interpolation of the list gives at the four
course frequencies compared with the exact model value.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from tissue import gabriel_props

COURSE = (403e6, 915e6, 2450e6, 5800e6)
TISSUE_OF = {"muscle": "muscle", "skin": "skin", "fat": "fat"}


def load(path):
    rows = []
    for r in csv.reader(open(path)):
        try:
            rows.append([float(x) for x in r[:3]])
        except ValueError:
            continue
    a = np.array(rows)
    if a.ndim != 2 or a.shape[1] != 3:
        raise SystemExit(f"{path}: expected at least 3 numeric columns")
    return a


def compare(paths):
    """paths: {tissue: file}. Returns {tissue: dict of statistics}."""
    out = {}
    for t, p in paths.items():
        a = load(p)
        f, er, sg = a[:, 0], a[:, 1], a[:, 2]
        mine = np.array([gabriel_props(TISSUE_OF[t], x) for x in f])
        de = np.abs(mine[:, 0] - er) / er
        ds = np.abs(mine[:, 1] - sg) / sg
        interp = []
        for fc in COURSE:
            e_l, s_l = np.interp(fc, f, er), np.interp(fc, f, sg)
            e_x, s_x = gabriel_props(TISSUE_OF[t], fc)
            interp.append((fc, e_l, e_x, s_l, s_x))
        out[t] = {"n": len(f), "fmin": f.min(), "fmax": f.max(), "dev_eps": de.max(), "dev_sigma": ds.max(),
                  "interp": interp, "data": a}
    return out


def write_3col(stats, out_dir):
    """Feko 'Import points' format: frequency, eps_r, sigma; full precision, no header."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for t, s in stats.items():
        with open(out_dir / f"{t}_course_3col.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            for f, er, sg in s["data"]:
                w.writerow([f"{f:.6e}", f"{er:.4f}", f"{sg:.6f}"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for t in TISSUE_OF:
        ap.add_argument(f"--{t}")
    ap.add_argument("--write-3col", metavar="DIR")
    args = ap.parse_args()
    paths = {t: getattr(args, t) for t in TISSUE_OF if getattr(args, t)}
    if not paths:
        raise SystemExit("give at least one of --muscle --skin --fat")
    stats = compare(paths)
    for t, s in stats.items():
        print(f"{t:7s} {s['n']} points {s['fmin']:.0e}-{s['fmax']:.0e} Hz: largest deviation "
              f"eps_r {100*s['dev_eps']:.4f} %, sigma {100*s['dev_sigma']:.4f} %")
    print("\nFeko linear interpolation of the list vs exact model at the course frequencies")
    print(f"{'tissue':7s}{'f MHz':>7}{'eps list':>10}{'eps exact':>10}{'sig list':>10}{'sig exact':>10}")
    for t, s in stats.items():
        for fc, e_l, e_x, s_l, s_x in s["interp"]:
            print(f"{t:7s}{fc/1e6:>7g}{e_l:>10.3f}{e_x:>10.3f}{s_l:>10.4f}{s_x:>10.4f}")
    if args.write_3col:
        write_3col(stats, args.write_3col)
        print(f"\nWrote 3-column files to {args.write_3col}")


if __name__ == "__main__":
    main()
