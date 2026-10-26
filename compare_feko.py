"""Week 5: compare a Feko near-field line with the script.

Usage:
  python compare_feko.py inputs/stack_2450_feko_slab.json feko_line.csv \
      [--zcol 0 --ecol 1 --zunit m] [--gamma 0.522] [--pabs skin=3.17,fat=0.55,muscle=3.55]

feko_line.csv  any text/CSV export of |E| along the line through the cell centre
               (POSTFEKO: right-click the 2D trace -> Export data). Non-numeric header
               lines are skipped; pick the z and |E| columns with --zcol/--ecol.
               Feko |E| must be the peak value for the same incident field as the JSON.
--gamma        |Gamma| read from Feko (reflection coefficient of the PBC solve)
--pabs         power loss per medium from Feko, in W/m^2 of cell area
               (Feko gives W per cell: divide by the 10 x 10 mm cell area, 1e-4 m^2)

Prints the checks of the Validation sheet: |Gamma| and per-layer absorption within 3 %,
muscle delta within 2 %, |Gamma| from the standing-wave ratio in air, and |E| error per layer.
Saves <feko stem>_vs_script.png next to the Feko file.
"""
import argparse
import re
from pathlib import Path

import numpy as np

from multilayer import fit_delta, solve
from run import load_case


def read_columns(path, zcol, ecol):
    z, e = [], []
    for line in Path(path).read_text().splitlines():
        parts = [p for p in re.split(r"[,\s;]+", line.strip()) if p]
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        if len(vals) > max(zcol, ecol):
            z.append(vals[zcol]); e.append(vals[ecol])
    if not z:
        raise SystemExit(f"no numeric rows found in {path}")
    z, e = np.array(z), np.array(e)
    order = np.argsort(-z)
    return z[order], e[order]


def gamma_from_swr(z, absE, f_hz):
    """|Gamma| = (S - 1)/(S + 1) from |E| in air (z > 0), using at least half a wavelength."""
    sel = z > 1e-3
    if np.ptp(z[sel]) < 0.5 * 3e8 / f_hz:
        return None
    s = absE[sel].max() / absE[sel].min()
    return (s - 1) / (s + 1)


def pct(a, b):
    return 100 * (a - b) / b


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case"); ap.add_argument("feko")
    ap.add_argument("--zcol", type=int, default=0); ap.add_argument("--ecol", type=int, default=1)
    ap.add_argument("--zunit", choices=["m", "mm"], default="m")
    ap.add_argument("--gamma", type=float)
    ap.add_argument("--pabs", help="e.g. skin=3.17,fat=0.55,muscle=3.55 (W/m^2)")
    args = ap.parse_args()

    cfg, f_hz, layers = load_case(args.case)
    res = solve(layers, f_hz, cfg.get("S_inc_Wpm2", 10.0))
    z, eF = read_columns(args.feko, args.zcol, args.ecol)
    if args.zunit == "mm":
        z = z * 1e-3
    eS = np.abs(res.fields(z)[0])

    print(f"Script: |Gamma| = {abs(res.Gamma):.5f}   Feko line: {len(z)} points, z = {z.max()*1e3:.1f} ... {z.min()*1e3:.1f} mm")
    if args.gamma is not None:
        e = pct(args.gamma, abs(res.Gamma))
        print(f"|Gamma| Feko {args.gamma:.5f} vs script: {e:+.2f} %  -> {'PASS' if abs(e) < 3 else 'FAIL'} (3 %)")
    g_swr = gamma_from_swr(z, eF, f_hz)
    if g_swr is not None:
        print(f"|Gamma| from SWR in air (Feko line): {g_swr:.4f} ({pct(g_swr, abs(res.Gamma)):+.2f} %)")

    print("\n|E| difference Feko vs script per layer (interface points skipped):")
    edges = [(0.0, np.inf, "air")] + [(-l.x_bottom, -l.x_top, l.name) for l in res.layers]
    for lo, hi, name in edges:
        sel = (z > lo + 0.25e-3) & (z < hi - 0.25e-3) if name != "air" else z > 0.25e-3
        if sel.sum() == 0:
            continue
        rel = np.abs(eF[sel] - eS[sel]) / eS[sel]
        print(f"  {name:<8} mean {100*rel.mean():5.2f} %   max {100*rel.max():5.2f} %   ({sel.sum()} pts)")

    muscle = res.layers[-1]
    zb = -(muscle.x_bottom if muscle.d is not None else muscle.x_top + 0.06)
    print(f"\nMuscle delta (analytic {muscle.delta*1e3:.2f} mm, pass: 2 %):")
    fits = [("Week 1 window (5 mm off each face)", None)]
    if muscle.d is not None:
        fits.append(("back-face-safe (skip 1.5 delta at bottom)", 1.5 * muscle.delta))
    for label, sb in fits:
        try:
            dS = fit_delta(z, eS, -muscle.x_top, zb, skip_bottom=sb)
            dF = fit_delta(z, eF, -muscle.x_top, zb, skip_bottom=sb)
            e = pct(dF, muscle.delta)
            print(f"  {label}: Feko {dF*1e3:.2f} mm ({e:+.2f} %, {'PASS' if abs(e) < 2 else 'FAIL'})"
                  f"   same fit on exact script field: {pct(dS, muscle.delta):+.2f} %")
        except ValueError as err:
            print(f"  {label}: {err}")

    if args.pabs:
        print("\nAbsorbed power per layer (W/m^2):")
        feko = dict(kv.split("=") for kv in args.pabs.split(","))
        for l in res.layers:
            if l.name in feko:
                v = float(feko[l.name]); e = pct(v, l.P_abs)
                print(f"  {l.name:<8} Feko {v:.4f}  script {l.P_abs:.4f}  {e:+.2f} %  -> {'PASS' if abs(e) < 3 else 'FAIL'}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    zz = np.linspace(z.max(), z.min(), 3000)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(zz * 1e3, np.abs(res.fields(zz)[0]), label="script (exact)")
    ax.plot(z * 1e3, eF, ".", ms=3, label="Feko")
    for l in res.layers:
        ax.axvline(-l.x_top * 1e3, color="0.6", lw=0.8, ls="--")
    ax.invert_xaxis(); ax.set_xlabel("z (mm)"); ax.set_ylabel("|E| (V/m, peak)")
    ax.set_title(f"{f_hz/1e6:g} MHz: Feko vs transmission-line script"); ax.grid(alpha=0.3); ax.legend()
    png = Path(args.feko).with_name(Path(args.feko).stem + "_vs_script.png")
    fig.tight_layout(); fig.savefig(png, dpi=150); plt.close(fig)
    print(f"\nSaved {png}")


if __name__ == "__main__":
    main()
