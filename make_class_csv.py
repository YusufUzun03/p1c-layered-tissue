"""Feko .out (near-field line on the cell/slab axis) -> class CSV schema.

Usage:
  python make_class_csv.py feko/C2.out --case plane_wave_feko --team T3 --out results/feko_plane_wave_2450.csv
  python make_class_csv.py feko/E3.out --case dipole_feko --team T3 --out results/feko_dipole_2450.csv

Columns: team, f_MHz, case, z_mm, layer, absE_Vpm, pabs_Wpm3, with absE the peak |E| of the total field
(Ex, Ey, Ez combined) and pabs = 1/2 sigma |E|^2 using the sigma that Feko reports for the medium of
the field point (the medium index column of the .out file). z_mm is the position on the line, air z > 0.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from parse_feko_out import parse


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--case", required=True)
    ap.add_argument("--team", default="T3")
    ap.add_argument("--out-file", "--out", dest="dest", required=True)
    args = ap.parse_args()

    d = parse(Path(args.out).read_text(errors="replace"))
    media = {m["index"]: m for m in d["media"]}
    rows = [r for r in d["near"] if abs(r[1]) < 1e-12 and abs(r[2]) < 1e-12]
    rows.sort(key=lambda r: -r[3])
    Path(args.dest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.dest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["team", "f_MHz", "case", "z_mm", "layer", "absE_Vpm", "pabs_Wpm3"])
        for med, x, y, z, ex, ey, ez in rows:
            e = float(np.sqrt(ex ** 2 + ey ** 2 + ez ** 2))
            m = media.get(med)
            name = "air" if med == 0 or m is None else m["name"].replace("_gabriel", "").replace("Muscle", "muscle")
            sigma = 0.0 if med == 0 or m is None else m["sigma"]
            w.writerow([args.team, f"{d['f_Hz']/1e6:g}", args.case, f"{z*1e3:.3f}", name, f"{e:.4f}",
                        f"{0.5*sigma*e**2:.5g}"])
    print(f"Wrote {len(rows)} points to {args.dest}")


if __name__ == "__main__":
    main()
