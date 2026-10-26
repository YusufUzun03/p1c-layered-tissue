"""Input file -> results. Usage:  python run.py inputs/stack_2450.json [--out results] [--no-plot]

Writes to the output folder:
  <stem>_summary.json   Gamma, power entering, absorbed power per layer, checks
  <stem>_profile.csv    class CSV schema: team, f_MHz, case, z_mm, layer, absE_Vpm, pabs_Wpm3
  <stem>_profile.png    |E|(z) and 1/2 sigma |E|^2 (z)
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from multilayer import make_layer, solve


def load_case(path):
    with open(path) as fh:
        cfg = json.load(fh)
    f_hz = cfg["f_MHz"] * 1e6
    layers = [make_layer(spec, f_hz) for spec in cfg["layers"]]
    return cfg, f_hz, layers


def profile_grid(res, cfg):
    """z grid [m]: from air_mm above the surface to the bottom of the stack, dz_mm step."""
    dz = cfg.get("dz_mm", 0.2) * 1e-3
    last = res.layers[-1]
    if last.d is None:
        depth = last.x_top + cfg.get("halfspace_depth_mm", 60.0) * 1e-3
    else:
        depth = last.x_bottom
    z_top = cfg.get("air_mm", 50.0) * 1e-3
    n = int(round((z_top + depth) / dz))
    return z_top - dz * np.arange(n + 1)


def print_summary(s):
    print(f"f = {s['f_MHz']:.0f} MHz   S_inc = {s['S_inc_Wpm2']} W/m^2   |E0| = {s['E0_Vpm_peak']:.2f} V/m (peak)")
    print(f"Gamma = {s['Gamma_re']:+.5f} {s['Gamma_im']:+.5f}j   |Gamma| = {s['abs_Gamma']:.5f}")
    print(f"Power entering tissue: {s['P_in_Wpm2']:.4f} W/m^2 ({100*s['P_in_Wpm2']/s['S_inc_Wpm2']:.2f} % of incident)")
    print(f"{'layer':<10}{'d (mm)':>9}{'eps_r':>9}{'sigma':>9}{'delta mm':>10}{'P_abs W/m2':>12}{'% inc':>8}{'% in':>8}")
    for l in s["layers"]:
        d = "half" if l["d_mm"] is None else f"{l['d_mm']:.2f}"
        print(f"{l['layer']:<10}{d:>9}{l['eps_r']:>9.2f}{l['sigma_Spm']:>9.3f}{l['delta_mm']:>10.1f}"
              f"{l['P_abs_Wpm2']:>12.4f}{100*l['frac_of_incident']:>8.2f}{100*l['frac_of_entering']:>8.2f}")
    print(f"Transmitted out of the bottom: {s['P_trans_Wpm2']:.4g} W/m^2")
    print(f"Power conservation error: {s['conservation_error']:.2e}   "
          f"Poynting vs integral(1/2 sigma|E|^2): {s['integral_error']:.2e}")
    for note in s["notes"]:
        print("NOTE:", note)


def write_profile(res, cfg, z, path):
    E, _, idx = res.fields(z)
    pabs = res.loss_density(z)
    names = [res.layers[i].name if 0 <= i < len(res.layers) else "air" for i in idx]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["team", "f_MHz", "case", "z_mm", "layer", "absE_Vpm", "pabs_Wpm3"])
        for zi, Ei, ni, pi in zip(z, E, names, pabs):
            w.writerow([cfg.get("team", ""), f"{res.f/1e6:g}", cfg.get("case", "plane_wave"),
                        f"{zi*1e3:.3f}", ni, f"{abs(Ei):.6g}", f"{pi:.6g}"])


def plot_profile(res, z, path, title=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    E, _, _ = res.fields(z)
    pabs = res.loss_density(z)
    zmm = z * 1e3
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax1.plot(zmm, np.abs(E), color="C0")
    ax1.set_ylabel("|E| (V/m, peak)")
    ax2.plot(zmm, pabs, color="C3")
    ax2.set_ylabel(r"$\frac{1}{2}\sigma|E|^2$ (W/m$^3$)")
    ax2.set_xlabel("z (mm)  — air z > 0, tissue z < 0")
    for ax in (ax1, ax2):
        for l in res.layers:
            ax.axvline(-l.x_top * 1e3, color="0.6", lw=0.8, ls="--")
        ax.grid(alpha=0.3)
    for l in res.layers:
        ax1.text(-l.x_top * 1e3, ax1.get_ylim()[1], f" {l.name}", va="top", ha="right",
                 rotation=90, fontsize=8, color="0.4")
    ax1.invert_xaxis()  # depth increases to the right
    ax1.set_title(title or f"{res.f/1e6:g} MHz, |Γ| = {abs(res.Gamma):.3f}")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="JSON case file")
    ap.add_argument("--out", default="results", help="output folder")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    cfg, f_hz, layers = load_case(args.input)
    res = solve(layers, f_hz, cfg.get("S_inc_Wpm2", 10.0))
    summary = res.summary()
    print_summary(summary)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(args.input).stem
    with open(out / f"{stem}_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    z = profile_grid(res, cfg)
    write_profile(res, cfg, z, out / f"{stem}_profile.csv")
    if not args.no_plot:
        plot_profile(res, z, out / f"{stem}_profile.png")
    print(f"\nWrote {out/stem}_summary.json, _profile.csv" + ("" if args.no_plot else ", _profile.png"))


if __name__ == "__main__":
    main()
