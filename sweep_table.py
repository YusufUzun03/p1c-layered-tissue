"""Zin, S11, accepted power and losses per frequency from a Feko .out file with several frequencies.

Usage:  python sweep_table.py feko/E4_sweep_partial.out [more.out ...] --csv results/dipole_sweep.csv
Blocks without an impedance (a run that aborted) are skipped.
"""
import argparse
import csv
import re
from pathlib import Path

F = r"[-+]?\d+\.?\d*(?:[Ee][-+]?\d+)?"


def blocks(text):
    parts = re.split(r"EXCITATION BY VOLTAGE SOURCE AT A SEGMENT", text)[1:]
    for p in parts:
        f = re.search(r"FREQ =\s*(" + F + ")", p)
        z = re.search(r"Impedance in Ohm\s+(" + F + r")\s+(" + F + ")", p)
        s = re.search(r"\n S\s+1\s+1\s+" + F + r"\s+" + F + r"\s+" + F + r"\s+(" + F + ")", p)
        pw = re.search(r"Power in Watt:\s*(" + F + ")", p)
        eff = re.search(r"Efficiency of the antenna:\s*(" + F + ")", p)
        if not (f and z and s and pw):
            continue
        loss = {}
        for n in ("skin_gabriel", "fat_gabriel", "Muscle_gabriel"):
            m = re.search(r"^\s*" + n + r"\s+.*?(" + F + r")\s*$", p, re.M)
            if m:
                loss[n.replace("_gabriel", "").replace("Muscle", "muscle")] = float(m.group(1))
        yield {"f_GHz": float(f.group(1)) / 1e9, "R": float(z.group(1)), "X": float(z.group(2)),
               "S11_dB": float(s.group(1)), "Pacc_W": float(pw.group(1)),
               "eff_pct": float(eff.group(1)) if eff else float("nan"), **loss}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", nargs="+")
    ap.add_argument("--csv")
    a = ap.parse_args()
    rows = []
    for o in a.out:
        rows += list(blocks(Path(o).read_text(errors="replace")))
    rows.sort(key=lambda r: r["f_GHz"])
    print(f"{'f GHz':>7}{'R':>8}{'X':>8}{'S11 dB':>9}{'Pacc mW':>9}{'eff %':>7}")
    for r in rows:
        print(f"{r['f_GHz']:>7.3f}{r['R']:>8.2f}{r['X']:>8.2f}{r['S11_dB']:>9.2f}{1e3*r['Pacc_W']:>9.3f}{r['eff_pct']:>7.1f}")
    if a.csv:
        Path(a.csv).parent.mkdir(parents=True, exist_ok=True)
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
