"""Week 6: dipole over the layered slab vs the plane wave.

Usage:
  python compare_dipole.py inputs/stack_2450_feko_slab.json --axis dipole_axis.csv \
      --pacc 0.1 --prad 0.012 --loss skin=0.030,fat=0.004,muscle=0.050 \
      [--fres-free 2.45 --fres-tissue 2.33] [--zcol 0 --ecol 1 --zunit m]

Inputs (all from POSTFEKO / the .out file of the finite-slab dipole run):
  --axis    |E| (peak) along the line below the feed, z = 0 at the skin surface, z < 0 in tissue
  --pacc    accepted power P_acc [W];  --prad  radiated power P_rad [W] (far-field integral)
  --loss    power loss per medium [W], layers named as in the JSON case
  --fres-*  resonance of the free-space dipole and of the dipole over tissue [GHz] (optional)

Prints the Week 3 / Validation-sheet checks and the comparison with the plane wave:
  * power balance  P_acc = P_rad + P_abs  (error < 5 %) and e_rad = P_rad / P_acc
  * fraction of the absorbed power per layer, dipole vs plane wave
  * on-axis depth profile of 1/2 sigma |E|^2, each case normalised by its own depth integral
    (so both describe where along the axis the power goes, independent of the source level)
    and the depth that holds 50 % / 90 % of it
Caveat: the on-axis profile of the dipole is a line through the maximum, not the total absorbed
power distribution. The per-layer fractions come from the power loss per medium, which is total.
"""
import argparse
from pathlib import Path

import numpy as np

from compare_feko import read_columns
from multilayer import solve
from run import load_case


def sigma_profile(res, z):
    sig = np.zeros(z.shape)
    for l in res.layers:
        bottom = l.x_bottom
        sel = (-z >= l.x_top) & (-z < bottom)
        sig[sel] = l.sigma
    return sig


def cumulative_depths(depth, p, fractions=(0.5, 0.9)):
    c = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(depth))])
    c = c / c[-1]
    return [float(np.interp(f, c, depth)) for f in fractions]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--axis", required=True)
    ap.add_argument("--pacc", type=float, required=True)
    ap.add_argument("--prad", type=float, required=True)
    ap.add_argument("--loss", required=True, help="skin=..,fat=..,muscle=.. in W")
    ap.add_argument("--fres-free", type=float)
    ap.add_argument("--fres-tissue", type=float)
    ap.add_argument("--zcol", type=int, default=0)
    ap.add_argument("--ecol", type=int, default=1)
    ap.add_argument("--zunit", choices=["m", "mm"], default="m")
    args = ap.parse_args()

    cfg, f_hz, layers = load_case(args.case)
    res = solve(layers, f_hz, cfg.get("S_inc_Wpm2", 10.0))

    loss = {k: float(v) for k, v in (kv.split("=") for kv in args.loss.split(","))}
    p_abs = sum(loss.values())
    err = abs(args.pacc - args.prad - p_abs) / args.pacc
    print(f"Power balance: P_acc = {args.pacc:.4g} W, P_rad = {args.prad:.4g} W, P_abs = {p_abs:.4g} W")
    print(f"  P_acc - P_rad - P_abs = {args.pacc - args.prad - p_abs:+.3g} W -> error {100*err:.2f} %"
          f"  -> {'PASS' if err < 0.05 else 'FAIL'} (< 5 %)")
    print(f"  radiation efficiency e_rad = P_rad / P_acc = {100*args.prad/args.pacc:.1f} %"
          f"   absorbed fraction = {100*p_abs/args.pacc:.1f} %")

    if args.fres_free and args.fres_tissue:
        shift = (args.fres_tissue - args.fres_free) / args.fres_free
        print(f"\nDetuning: f_res {args.fres_free:.3f} GHz (free space) -> {args.fres_tissue:.3f} GHz over tissue"
              f" ({100*shift:+.1f} %)")

    print("\nShare of the absorbed power per layer")
    print(f"  {'layer':<8}{'dipole':>9}{'plane wave':>12}{'difference':>12}")
    for l in res.layers:
        d = loss.get(l.name)
        pw = l.P_abs / sum(x.P_abs for x in res.layers)
        if d is None:
            continue
        print(f"  {l.name:<8}{100*d/p_abs:>8.1f}%{100*pw:>11.1f}%{100*(d/p_abs-pw):>+11.1f}")
    missing = [l.name for l in res.layers if l.name not in loss]
    if missing:
        print(f"  (no Feko loss given for: {', '.join(missing)})")

    z, e = read_columns(args.axis, args.zcol, args.ecol)
    if args.zunit == "mm":
        z = z * 1e-3
    sel = z <= 0
    z, e = z[sel], e[sel]
    order = np.argsort(-z)
    z, e = z[order], e[order]
    depth = -z
    p_dip = 0.5 * sigma_profile(res, z) * e ** 2
    p_pw = res.loss_density(z)
    n_dip = p_dip / np.trapezoid(p_dip, depth)
    n_pw = p_pw / np.trapezoid(p_pw, depth)
    d_dip = cumulative_depths(depth, p_dip)
    d_pw = cumulative_depths(depth, p_pw)
    print("\nOn-axis absorbed-power profile (each case normalised by its own integral over depth)")
    print(f"  depth holding 50 % / 90 %:  dipole {d_dip[0]*1e3:.1f} / {d_dip[1]*1e3:.1f} mm"
          f"   plane wave {d_pw[0]*1e3:.1f} / {d_pw[1]*1e3:.1f} mm")
    print(f"  surface value (first point in tissue): dipole {n_dip[0]*1e-3:.1f} /mm, plane wave {n_pw[0]*1e-3:.1f} /mm")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(depth * 1e3, n_pw * 1e-3, label="plane wave (script)")
    ax.plot(depth * 1e3, n_dip * 1e-3, label="dipole, on-axis (Feko)")
    for l in res.layers:
        ax.axvline(l.x_top * 1e3, color="0.7", lw=0.8, ls="--")
    ax.set_xlabel("depth below the skin surface (mm)")
    ax.set_ylabel("½σ|E|², normalised (1/mm)")
    ax.set_title(f"{f_hz/1e6:g} MHz: where the power goes, dipole vs plane wave")
    ax.grid(alpha=0.3); ax.legend()
    png = Path(args.axis).with_name(Path(args.axis).stem + "_vs_planewave.png")
    fig.tight_layout(); fig.savefig(png, dpi=150); plt.close(fig)
    print(f"\nSaved {png}")


if __name__ == "__main__":
    main()
