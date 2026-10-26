"""Checkpoint 1 / Week 5: turn Feko runs into model-log rows and convergence evidence.

Usage:
  python model_log.py feko/model_b_log.csv                     # Model B: muscle slab T = 60 mm, 2.45 GHz
  python model_log.py feko/model_b_log.csv --T-mm 90           # truncation check run set
  python model_log.py feko/stack_log.csv --case inputs/stack_2450_feko_slab.json --key abs_E

Log CSV (one row per Feko solve, see feko/model_b_log_template.csv):
  run, date, change, solver, mesh_mm, unknowns, mem_GB, time_min, gamma, line_file, notes
  - gamma      |Gamma| from Feko (optional if line_file is given: then taken from the SWR in air)
  - line_file  POSTFEKO export of |E| along the centre line (z in m, |E| peak), path relative
               to the log file; used for delta and |E| at the probe point. Optional.

Writes next to the log:
  <log>_workbook.csv     rows in the column order of the workbook's 'Model log' sheet
  <log>_convergence.png  each quantity vs mesh size, with the script reference
Checks printed: convergence between the two finest meshes (< 2 % fields, Week 4),
agreement with the script (|Gamma| and delta < 2 %, Week 1), laptop budget (12 GB, 15 min).
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from compare_feko import gamma_from_swr, read_columns
from multilayer import fit_delta, make_layer, solve
from run import load_case
from tissue import ETA0

UNITS = {"gamma": "-", "delta_mm": "mm", "abs_E": "V/m"}
LABELS = {"gamma": "|Gamma|", "delta_mm": "delta (muscle)", "abs_E": "|E| at probe"}


def reference(args):
    if args.case:
        cfg, f_hz, layers = load_case(args.case)
        S = cfg.get("S_inc_Wpm2", 10.0)
    else:
        f_hz, S = args.f_MHz * 1e6, args.S_inc
        layers = [make_layer({"name": "muscle", "tissue": "muscle", "d_mm": args.T_mm}, f_hz)]
    res = solve(layers, f_hz, S)
    m = res.layers[-1]
    z_probe = -(m.x_top + args.probe_mm * 1e-3)
    return res, m, z_probe


def quantities(row, base, res, m, z_probe):
    q = {"gamma": float(row["gamma"]) if row.get("gamma") else None, "delta_mm": None, "abs_E": None}
    if row.get("line_file"):
        z, e = read_columns(base / row["line_file"], 0, 1)
        if q["gamma"] is None:
            q["gamma"] = gamma_from_swr(z, e, res.f)
        zb = -(m.x_bottom if m.d is not None else m.x_top + 0.06)
        sb = 1.5 * m.delta if m.d is not None else None
        try:
            q["delta_mm"] = fit_delta(z, e, -m.x_top, zb, skip_bottom=sb) * 1e3
        except ValueError:
            pass
        if z.min() <= z_probe <= z.max():
            q["abs_E"] = float(np.interp(z_probe, z[::-1], e[::-1]))
    return q


def rel(a, b):
    return None if a is None or b is None else (a - b) / b


def fmt_pct(x):
    return "   -   " if x is None else f"{100*x:+6.2f}%"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log")
    ap.add_argument("--case", help="JSON case instead of the built-in Model B slab")
    ap.add_argument("--f-MHz", type=float, default=2450.0)
    ap.add_argument("--T-mm", type=float, default=60.0, help="muscle slab thickness for Model B")
    ap.add_argument("--S-inc", type=float, default=10.0, help="incident W/m^2 used in Feko")
    ap.add_argument("--E0", type=float, help="incident |E| peak in V/m used in Feko (e.g. 1); overrides --S-inc")
    ap.add_argument("--probe-mm", type=float, default=10.0, help="probe depth below top of muscle")
    ap.add_argument("--key", choices=list(UNITS), default="gamma", help="key output for the workbook rows")
    args = ap.parse_args()
    if args.E0 is not None:
        args.S_inc = args.E0 ** 2 / (2 * ETA0)

    log = Path(args.log)
    with open(log, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r.get("run", "").strip()]
    if not rows:
        raise SystemExit("log has no runs yet")

    res, m, z_probe = reference(args)
    ref = {"gamma": abs(res.Gamma), "delta_mm": m.delta * 1e3, "abs_E": float(np.abs(res.fields(z_probe)[0][0]))}
    qs = [quantities(r, log.parent, res, m, z_probe) for r in rows]

    print(f"Reference (script): |Gamma| = {ref['gamma']:.5f}, delta = {ref['delta_mm']:.2f} mm, "
          f"|E| at {args.probe_mm:g} mm into muscle = {ref['abs_E']:.3f} V/m (S_inc = {res.S_inc} W/m^2)\n")
    print(f"{'run':<5}{'mesh':>7}{'unknowns':>10}{'GB':>6}{'min':>6}"
          f"{'|Gamma|':>9}{'chg':>9}{'vs ref':>9}{'delta':>8}{'chg':>9}{'vs ref':>9}{'|E|pr':>8}{'chg':>9}")
    prev = {k: None for k in UNITS}
    changes = []
    for r, q in zip(rows, qs):
        ch = {k: rel(q[k], prev[k]) for k in UNITS}
        changes.append(ch)
        g = "      - " if q["gamma"] is None else f"{q['gamma']:>8.4f} "
        d = "     - " if q["delta_mm"] is None else f"{q['delta_mm']:>7.2f} "
        e = "     - " if q["abs_E"] is None else f"{q['abs_E']:>7.3f} "
        print(f"{r['run']:<5}{r.get('mesh_mm',''):>7}{r.get('unknowns',''):>10}{r.get('mem_GB',''):>6}{r.get('time_min',''):>6}"
              f"{g}{fmt_pct(ch['gamma']):>9}{fmt_pct(rel(q['gamma'], ref['gamma'])):>9}"
              f"{d}{fmt_pct(ch['delta_mm']):>9}{fmt_pct(rel(q['delta_mm'], ref['delta_mm'])):>9}"
              f"{e}{fmt_pct(ch['abs_E']):>9}")
        prev = {k: q[k] if q[k] is not None else prev[k] for k in UNITS}

    print("\nChecks")
    last, ch_last = qs[-1], changes[-1]
    for k, tol in (("gamma", 0.02), ("delta_mm", 0.02), ("abs_E", 0.02)):
        if ch_last[k] is not None:
            print(f"  convergence {LABELS[k]:<16} last refinement {fmt_pct(ch_last[k])}  "
                  f"-> {'PASS' if abs(ch_last[k]) < tol else 'FAIL'} (< {100*tol:g} %)")
    for k in ("gamma", "delta_mm"):
        e = rel(last[k], ref[k])
        if e is not None:
            print(f"  vs script   {LABELS[k]:<16} finest mesh {fmt_pct(e)}  -> {'PASS' if abs(e) < 0.02 else 'FAIL'} (< 2 %)")
    mem = [float(r["mem_GB"]) for r in rows if r.get("mem_GB")]
    tim = [float(r["time_min"]) for r in rows if r.get("time_min")]
    if mem or tim:
        ok = (not mem or max(mem) <= 12) and (not tim or max(tim) <= 15)
        print(f"  laptop budget: peak {max(mem) if mem else '-'} GB, longest {max(tim) if tim else '-'} min"
              f"  -> {'PASS' if ok else 'FAIL'} (12 GB, 15 min)")
    if len(rows) < 3:
        print("  NOTE: Week 4 asks for at least 3 meshes (lambda/10, /15, /20)")

    out_csv = log.with_name(log.stem + "_workbook.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Run", "Date", "Change from previous run", "Solver", "Mesh (edge)", "Unknowns / tets",
                    "Peak memory (GB)", "Time (min)", "Key output", "Value", "Unit", "Rel. change vs previous", "Notes"])
        for r, q, ch in zip(rows, qs, changes):
            v = q[args.key]
            w.writerow([r["run"], r.get("date", ""), r.get("change", ""), r.get("solver", ""),
                        f"{r['mesh_mm']} mm" if r.get("mesh_mm") else "", r.get("unknowns", ""),
                        r.get("mem_GB", ""), r.get("time_min", ""), LABELS[args.key],
                        "" if v is None else f"{v:.5g}", UNITS[args.key],
                        "" if ch[args.key] is None else f"{100*ch[args.key]:.2f} %", r.get("notes", "")])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mesh = [float(r["mesh_mm"]) if r.get("mesh_mm") else np.nan for r in rows]
    keys = [k for k in UNITS if any(q[k] is not None for q in qs)]
    fig, axes = plt.subplots(1, len(keys), figsize=(4 * len(keys), 3.4), squeeze=False)
    for ax, k in zip(axes[0], keys):
        pts = [(x, q[k]) for x, q in zip(mesh, qs) if q[k] is not None and not np.isnan(x)]
        if pts:
            ax.plot(*zip(*pts), "o-", label="Feko")
        ax.axhline(ref[k], ls="--", color="0.5", label="script")
        ax.set_xlabel("mesh edge (mm)"); ax.set_title(f"{LABELS[k]} ({UNITS[k]})"); ax.invert_xaxis()
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout()
    png = log.with_name(log.stem + "_convergence.png")
    fig.savefig(png, dpi=150); plt.close(fig)
    print(f"\nWrote {out_csv} and {png}")


if __name__ == "__main__":
    main()
