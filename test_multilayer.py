"""Validation tests for Checkpoint 1.  Run:  pytest -q   or   python test_multilayer.py  (prints a report)

Test 1  single air/muscle interface vs Fresnel, 4 significant digits
Test 2  lossless lambda/4 matching layer, |Gamma| < 1e-6
Plus    power conservation, Gabriel media vs the Week 1 table (1 %),
        penetration depth from the field decay (2 %), field continuity.
"""
import numpy as np

from multilayer import Layer, fit_delta, make_layer, solve
from tissue import complex_eps, gabriel_eps, gabriel_props, wave_params

F = 2.45e9
SPEC_STACK = [{"tissue": "skin", "d_mm": 1.5}, {"tissue": "fat", "d_mm": 8}, {"tissue": "muscle", "d_mm": None}]

# Week 1 slide 17 (eps_r, sigma) at 403 / 915 / 2450 / 5800 MHz
WEEK1_TABLE = {
    "muscle": [(57.1, 0.797), (55.0, 0.948), (52.7, 1.739), (48.5, 4.962)],
    "skin": [(46.7, 0.689), (41.3, 0.872), (38.0, 1.464), (35.1, 3.717)],
    "fat": [(5.6, 0.041), (5.5, 0.051), (5.3, 0.105), (5.0, 0.293)],
}
FREQS = [403e6, 915e6, 2450e6, 5800e6]


def sig_digits_equal(a, b, n=4):
    return f"{a:.{n - 1}e}" == f"{b:.{n - 1}e}"


def stack(specs, f=F):
    return [make_layer(s, f) for s in specs]


def run_fresnel():
    """Air/muscle half-space vs Fresnel Gamma = (1 - n)/(1 + n)."""
    out = []
    for f in FREQS:
        res = solve(stack([{"tissue": "muscle", "d_mm": None}], f), f)
        n = np.sqrt(gabriel_eps("muscle", f))
        g_ref = (1 - n) / (1 + n)
        ok = sig_digits_equal(abs(res.Gamma), abs(g_ref)) and abs(res.Gamma - g_ref) < 1e-9
        out.append((f, abs(res.Gamma), abs(g_ref), 1 - abs(res.Gamma) ** 2, ok))
    return out


def run_quarter_wave():
    """Lossless lambda/4 layer with eta1 = sqrt(eta0 eta2) on a lossless half-space: Gamma = 0."""
    eps2 = 9.0
    eps1 = np.sqrt(eps2)                    # n1 = sqrt(n0 n2)
    lam1 = wave_params(complex_eps(eps1, 0, F), F)[3]
    layers = [Layer("match", lam1 / 4, complex_eps(eps1, 0, F), F),
              Layer("load", None, complex_eps(eps2, 0, F), F)]
    res = solve(layers, F)
    # off-design check: at 2f the layer is lambda/2 and Gamma equals the bare interface
    layers2 = [Layer("match", lam1 / 4, complex_eps(eps1, 0, 2 * F), 2 * F),
               Layer("load", None, complex_eps(eps2, 0, 2 * F), 2 * F)]
    res2 = solve(layers2, 2 * F)
    bare = abs((1 - np.sqrt(eps2)) / (1 + np.sqrt(eps2)))
    return abs(res.Gamma), abs(res2.Gamma), bare


def test_fresnel_4_digits():
    for f, g, g_ref, _, ok in run_fresnel():
        assert ok, f"{f/1e6} MHz: {g} vs {g_ref}"


def test_fresnel_matches_week1_slide():
    (_, _, _, _, _), (_, _, _, _, _), (f, g, _, t, _), _ = run_fresnel()
    assert round(g, 3) == 0.762 and round(100 * t, 1) == 41.9


def test_quarter_wave_match():
    g, g_half, bare = run_quarter_wave()
    assert g < 1e-6
    assert abs(g_half - bare) < 1e-9


def test_power_conservation():
    for f in FREQS:
        for specs in (SPEC_STACK, SPEC_STACK[:2] + [{"tissue": "muscle", "d_mm": 60}]):
            res = solve(stack(specs, f), f)
            assert res.conservation_error < 1e-6
            assert res.integral_error < 1e-4
            # entering power from Gamma equals Poynting flux just inside the skin
            assert abs(res.poynting(-1e-9)[0] - res.P_in) / res.P_in < 1e-6


def test_field_continuity():
    res = solve(stack(SPEC_STACK), F)
    for l in res.layers[1:]:
        z = -l.x_top
        E_above, H_above, _ = res.fields(z + 1e-12)
        E_below, H_below, _ = res.fields(z - 1e-12)
        assert abs(E_above[0] - E_below[0]) < 1e-6 * abs(E_below[0])
        assert abs(H_above[0] - H_below[0]) < 1e-6 * abs(H_below[0])


def test_gabriel_vs_week1_table():
    for tissue, rows in WEEK1_TABLE.items():
        for f, (er, sig) in zip(FREQS, rows):
            er_c, sig_c = gabriel_props(tissue, f)
            # slide values are rounded to 2-3 digits: 1 % tolerance
            assert abs(er_c - er) / er < 0.01, (tissue, f, er_c, er)
            assert abs(sig_c - sig) / sig < 0.01, (tissue, f, sig_c, sig)


def test_penetration_depth_fit():
    res = solve(stack([{"tissue": "muscle", "d_mm": None}]), F)
    z = -np.arange(0, 0.06, 0.2e-3)
    E, _, _ = res.fields(z)
    d_fit = fit_delta(z, np.abs(E), 0.0, -0.06)
    d_ref = res.layers[0].delta
    assert abs(d_fit - d_ref) / d_ref < 0.02
    assert abs(d_ref * 100 - 2.23) < 0.01  # Week 1: 2.23 cm


def test_penetration_depth_fit_finite_slab():
    """60 mm muscle slab backed by air: back-face ripple; skipping 1.5 delta at the bottom fixes it."""
    specs = SPEC_STACK[:2] + [{"tissue": "muscle", "d_mm": 60}]
    res = solve(stack(specs), F)
    m = res.layers[-1]
    z = -np.arange(0, m.x_bottom, 0.2e-3)
    absE = np.abs(res.fields(z)[0])
    d_safe = fit_delta(z, absE, -m.x_top, -m.x_bottom, skip_bottom=1.5 * m.delta)
    assert abs(d_safe - m.delta) / m.delta < 0.01


def test_dipole_reference():
    """Induced-EMF dipole: Week 2 target 73.1 + j42.5 ohm at 0.5 lambda, resonance a little below."""
    from dipole_reference import resonant_length, zin
    z = zin(0.5, 1e-3)
    assert abs(z.real - 73.1) / 73.1 < 0.01 and abs(z.imag - 42.5) < 0.5
    x0 = resonant_length(1e-3)
    assert 0.47 < x0 < 0.49
    assert 62 < zin(x0, 1e-3).real < 66  # Week 2: R about 65 ohm at resonance


def report():
    print("Checkpoint 1 validation report\n" + "=" * 31)
    print("\nTest 1: air/muscle interface vs Fresnel (pass: 4 significant digits)")
    print(f"{'f (MHz)':>8}{'script |G|':>14}{'Fresnel |G|':>14}{'1-|G|^2':>10}  pass")
    for f, g, g_ref, t, ok in run_fresnel():
        print(f"{f/1e6:>8.0f}{g:>14.6f}{g_ref:>14.6f}{t:>10.4f}  {'yes' if ok else 'NO'}")

    g, g_half, bare = run_quarter_wave()
    print("\nTest 2: lossless lambda/4 matching layer (eps1 = 3 on eps2 = 9), pass: |G| < 1e-6")
    print(f"  |G| at design frequency = {g:.2e}  -> {'yes' if g < 1e-6 else 'NO'}")
    print(f"  |G| at 2f (lambda/2 layer) = {g_half:.6f}, bare interface = {bare:.6f}")

    print("\nPower conservation: 1 - |G|^2 = sum P_abs/S_inc + P_trans/S_inc (spec stack)")
    print(f"{'f (MHz)':>8}{'|G|':>8}{'1-|G|^2':>10}{'sum abs':>10}{'trans':>10}{'error':>10}{'integral err':>14}")
    for f in FREQS:
        res = solve(stack(SPEC_STACK, f), f)
        tot = sum(l.P_abs for l in res.layers) / res.S_inc
        print(f"{f/1e6:>8.0f}{abs(res.Gamma):>8.4f}{1-abs(res.Gamma)**2:>10.5f}{tot:>10.5f}"
              f"{res.P_trans/res.S_inc:>10.2e}{res.conservation_error:>10.1e}{res.integral_error:>14.1e}")

    res = solve(stack([{"tissue": "muscle", "d_mm": None}]), F)
    z = -np.arange(0, 0.06, 0.2e-3)
    d_fit = fit_delta(z, np.abs(res.fields(z)[0]), 0.0, -0.06)
    print(f"\nMuscle delta at 2.45 GHz: analytic {res.layers[0].delta*1e3:.2f} mm, "
          f"fit of ln|E| {d_fit*1e3:.2f} mm")

    print("\nGabriel media vs Week 1 table (pass: 1 %)")
    worst = 0
    for tissue, rows in WEEK1_TABLE.items():
        for f, (er, sig) in zip(FREQS, rows):
            er_c, sig_c = gabriel_props(tissue, f)
            worst = max(worst, abs(er_c - er) / er, abs(sig_c - sig) / sig)
    print(f"  worst deviation {100*worst:.2f} %")


if __name__ == "__main__":
    report()


def test_parse_feko_out_fixture():
    from pathlib import Path
    from parse_feko_out import parse, derived
    d = parse((Path(__file__).parent / "tests" / "data" / "modelB_sample.out").read_text())
    x = derived(d)
    assert d["dof"] == 59628 and d["tets"] == 9984
    assert abs(d["R"][0] - 0.759806) < 1e-9
    assert abs(x["absorbed_fraction"] - 0.42469) < 1e-4
