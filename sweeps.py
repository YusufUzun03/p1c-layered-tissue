"""Parameter sweeps for the report.  Usage:  python sweeps.py [--out results]

  freq_sweep.csv / .png   Gamma and absorbed fraction per layer at 403, 915, 2450, 5800 MHz
                          (plus a dense 0.3-6 GHz curve for Output 1: Gamma vs frequency)
  fat_sweep.csv / .png    fat thickness 0-30 mm at 2.45 GHz: |Gamma|, absorbed power per layer
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from multilayer import make_layer, solve

COURSE_FREQS = [403e6, 915e6, 2450e6, 5800e6]


def spec_stack(fat_mm=8.0, skin_mm=1.5):
    return [{"name": "skin", "tissue": "skin", "d_mm": skin_mm},
            {"name": "fat", "tissue": "fat", "d_mm": fat_mm},
            {"name": "muscle", "tissue": "muscle", "d_mm": None}]


def run(specs, f, S_inc=10.0):
    return solve([make_layer(s, f) for s in specs], f, S_inc)


def freq_sweep(out):
    rows = []
    for f in COURSE_FREQS:
        r = run(spec_stack(), f)
        rows.append([f / 1e6, abs(r.Gamma), r.P_in / r.S_inc] + [l.P_abs / r.P_in for l in r.layers])
    with open(out / "freq_sweep.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["f_MHz", "abs_Gamma", "frac_entering", "skin_frac_of_in", "fat_frac_of_in", "muscle_frac_of_in"])
        w.writerows([[f"{v:.6g}" for v in row] for row in rows])

    fd = np.linspace(0.3e9, 6e9, 400)
    gd = [abs(run(spec_stack(), f).Gamma) for f in fd]
    g_muscle = [abs(run([{"tissue": "muscle", "d_mm": None}], f).Gamma) for f in fd]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(fd / 1e9, gd, label="skin 1.5 / fat 8 / muscle")
    ax1.plot(fd / 1e9, g_muscle, ls="--", color="0.5", label="muscle half-space")
    for row in rows:
        ax1.plot(row[0] / 1e3, row[1], "o", color="C0")
    ax1.set_xlabel("f (GHz)"); ax1.set_ylabel("|Γ|"); ax1.grid(alpha=0.3); ax1.legend()
    ax1.set_title("Output 1: reflection coefficient vs frequency")

    x = np.arange(len(rows))
    bottom = np.zeros(len(rows))
    for j, name in enumerate(["skin", "fat", "muscle"]):
        vals = np.array([row[3 + j] for row in rows]) * 100
        ax2.bar(x, vals, bottom=bottom, label=name)
        bottom += vals
    ax2.set_xticks(x, [f"{row[0]:g} MHz" for row in rows])
    ax2.set_ylabel("% of power entering tissue"); ax2.legend(loc="upper left", bbox_to_anchor=(1.01, 1)); ax2.set_title("Output 3: absorbed fraction per layer")
    fig.tight_layout(); fig.savefig(out / "freq_sweep.png", dpi=150); plt.close(fig)
    return rows


def fat_sweep(out, f=2.45e9, fats=np.arange(0, 30.01, 0.5)):
    rows = []
    for d in fats:
        specs = spec_stack(fat_mm=d) if d > 0 else [spec_stack()[0], spec_stack()[2]]
        r = run(specs, f)
        p = {l.name: l.P_abs for l in r.layers}
        rows.append([d, abs(r.Gamma), r.P_in, p["skin"], p.get("fat", 0.0), p["muscle"]])
    with open(out / "fat_sweep.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["fat_mm", "abs_Gamma", "P_in_Wpm2", "skin_Wpm2", "fat_Wpm2", "muscle_Wpm2"])
        w.writerows([[f"{v:.6g}" for v in row] for row in rows])

    a = np.array(rows)
    no_fat = a[0, 5]
    better = a[a[:, 5] > no_fat, 0]
    i_peak = a[:, 5].argmax()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(a[:, 0], a[:, 5], lw=2, label="muscle")
    ax.plot(a[:, 0], a[:, 3], label="skin")
    ax.plot(a[:, 0], a[:, 4], label="fat")
    ax.axhline(no_fat, ls="--", color="0.5", lw=1, label=f"muscle, no fat ({no_fat:.2f})")
    ax.axvline(8, ls=":", color="0.3", lw=1)
    ax.text(8.3, ax.get_ylim()[1] * 0.95, "spec 8 mm", va="top", fontsize=8)
    ax.set_xlabel("fat thickness (mm)"); ax.set_ylabel("absorbed power (W/m², 10 W/m² incident)")
    ax.set_title(f"Output 4: fat-thickness study at {f/1e9:g} GHz — muscle peak "
                 f"{a[i_peak,5]:.2f} W/m² at {a[i_peak,0]:g} mm")
    ax2 = ax.twinx(); ax2.plot(a[:, 0], a[:, 1], color="0.7", lw=1); ax2.set_ylabel("|Γ| (grey)")
    ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(out / "fat_sweep.png", dpi=150); plt.close(fig)
    return a, no_fat, better, i_peak


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results")
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)

    print("Frequency sweep (skin 1.5 / fat 8 / muscle half-space)")
    print(f"{'f MHz':>7}{'|G|':>8}{'enters':>9}{'skin':>8}{'fat':>8}{'muscle':>8}")
    for row in freq_sweep(out):
        print(f"{row[0]:>7.0f}{row[1]:>8.3f}{100*row[2]:>8.1f}%{100*row[3]:>7.1f}%{100*row[4]:>7.1f}%{100*row[5]:>7.1f}%")

    a, no_fat, better, i_peak = fat_sweep(out)
    print("\nFat sweep at 2.45 GHz")
    print(f"  |G| min {a[:,1].min():.3f} at {a[a[:,1].argmin(),0]:g} mm, max {a[:,1].max():.3f} at {a[a[:,1].argmax(),0]:g} mm")
    print(f"  muscle absorbs more than with no fat ({no_fat:.2f} W/m^2) for fat {better.min():g}-{better.max():g} mm;"
          f" peak {a[i_peak,5]:.2f} W/m^2 at {a[i_peak,0]:g} mm")
    print(f"\nWrote {out}/freq_sweep.csv/.png and {out}/fat_sweep.csv/.png")


if __name__ == "__main__":
    main()
