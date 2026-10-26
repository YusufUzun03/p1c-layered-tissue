"""Free-space thin-dipole reference for Week 6 (and the Week 2 Demo 1 self-test).

Induced-EMF theory with a sinusoidal current (Balanis, Antenna Theory, eqs. 8-60 to 8-61),
resistance and reactance referred to the input-terminal current:

    R_in = R_r / sin^2(kl/2),    X_in = X_m / sin^2(kl/2)

Usage:
  python dipole_reference.py                       # table, resonance length, plot at 2.45 GHz
  python dipole_reference.py --f-GHz 2.45 --a-over-lam 0.001 --out results

Targets for the Feko free-space dipole (Week 2): Z_in(0.5 lambda) = 73.1 + j42.5 ohm
(R within 3 %, X within 5 ohm), D = 2.15 dBi, resonance a little below 0.5 lambda.
Feko is a full-wave MoM solution with a gap and a finite wire, so it will not hit the
induced-EMF numbers exactly: use them as the reference the Week 2 pass criteria refer to.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import sici

from tissue import C0, ETA0

EULER = 0.5772156649015329


def zin(l_over_lam, a_over_lam=1e-3):
    """Input impedance (ohm) of a centre-fed dipole of total length l, radius a (both in wavelengths)."""
    x = np.asarray(l_over_lam, float)
    kl = 2 * np.pi * x
    Si, Ci = (lambda v: sici(v)[0]), (lambda v: sici(v)[1])
    arg = 4 * np.pi * a_over_lam ** 2 / x          # 2 k a^2 / l with lengths in wavelengths
    r_r = ETA0 / (2 * np.pi) * (
        EULER + np.log(kl) - Ci(kl)
        + 0.5 * np.sin(kl) * (Si(2 * kl) - 2 * Si(kl))
        + 0.5 * np.cos(kl) * (EULER + np.log(kl / 2) + Ci(2 * kl) - 2 * Ci(kl)))
    x_m = ETA0 / (4 * np.pi) * (
        2 * Si(kl) + np.cos(kl) * (2 * Si(kl) - Si(2 * kl))
        - np.sin(kl) * (2 * Ci(kl) - Ci(2 * kl) - Ci(arg)))
    s2 = np.sin(kl / 2) ** 2
    return r_r / s2 + 1j * x_m / s2


def resonant_length(a_over_lam=1e-3, lo=0.40, hi=0.50):
    """First length (in wavelengths) where X_in = 0."""
    return brentq(lambda v: zin(v, a_over_lam).imag, lo, hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--f-GHz", type=float, default=2.45)
    ap.add_argument("--a-over-lam", type=float, default=1e-3)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    lam = C0 / (args.f_GHz * 1e9)
    a = args.a_over_lam
    print(f"f = {args.f_GHz} GHz, lambda = {lam*1e3:.2f} mm, wire radius a = {a:g} lambda = {a*lam*1e3:.4f} mm\n")
    print(f"{'L/lambda':>9}{'L (mm)':>9}{'R (ohm)':>10}{'X (ohm)':>10}")
    for x in (0.40, 0.44, 0.46, 0.47, 0.48, 0.49, 0.50, 0.52):
        z = zin(x, a)
        print(f"{x:>9.2f}{x*lam*1e3:>9.2f}{z.real:>10.1f}{z.imag:>10.1f}")

    x0 = resonant_length(a)
    z0 = zin(x0, a)
    print(f"\nResonance (X = 0): L = {x0:.4f} lambda = {x0*lam*1e3:.2f} mm, R = {z0.real:.1f} ohm")
    z5 = zin(0.5, a)
    print(f"Half-wave (0.5 lambda = {0.5*lam*1e3:.2f} mm): Z = {z5.real:.1f} + j{z5.imag:.1f} ohm"
          f"  (Week 2 target 73.1 + j42.5)")
    gam = abs((z0 - 50) / (z0 + 50))
    print(f"With a 50 ohm port at resonance: |S11| = {20*np.log10(gam):.1f} dB")
    print("Directivity of a half-wave dipole: 1.64 = 2.15 dBi")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.linspace(0.40, 0.60, 400)
    z = zin(x, a)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, z.real, label="R_in")
    ax.plot(x, z.imag, label="X_in")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.axvline(x0, ls="--", color="0.5", lw=0.8)
    ax.annotate(f"X = 0 at {x0:.3f} λ\n({x0*lam*1e3:.1f} mm, R = {z0.real:.0f} Ω)", (x0, 0), (x0 + 0.02, -90),
                arrowprops=dict(arrowstyle="->", color="0.4"), fontsize=9)
    ax.plot([0.5], [z5.real], "o", color="C0"); ax.plot([0.5], [z5.imag], "o", color="C1")
    ax.set_xlabel("total length L / λ"); ax.set_ylabel("Ω"); ax.grid(alpha=0.3); ax.legend()
    ax.set_title(f"Thin dipole input impedance, a = {a:g} λ (induced-EMF reference)")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(out / "dipole_zin.png", dpi=150); plt.close(fig)
    print(f"\nWrote {out}/dipole_zin.png")


if __name__ == "__main__":
    main()
