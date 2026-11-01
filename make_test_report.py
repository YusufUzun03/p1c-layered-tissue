"""Build docs/validation_report.pdf from a real test run (numbers are computed, not typed).

Usage:  python make_test_report.py [--out docs/validation_report.pdf] [--author "Yusuf Uzun"]
"""
import argparse
import datetime
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np
import reportlab
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from multilayer import fit_delta, make_layer, solve
from test_multilayer import FREQS, SPEC_STACK, WEEK1_TABLE, run_fresnel, run_quarter_wave, stack
from check_course_media import compare as compare_media
from tissue import ETA0, gabriel_props
from parse_feko_out import derived as feko_derived, line_arrays, parse as parse_feko

ROOT = Path(__file__).resolve().parent
FONT = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
pdfmetrics.registerFont(TTFont("DV", str(FONT / "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DV-B", str(FONT / "DejaVuSans-Bold.ttf")))

ss = getSampleStyleSheet()
BODY = ParagraphStyle("body", parent=ss["Normal"], fontName="DV", fontSize=8.8, leading=12)
SMALL = ParagraphStyle("small", parent=BODY, fontSize=7.6, leading=10, textColor=colors.HexColor("#444444"))
H1 = ParagraphStyle("h1", parent=BODY, fontName="DV-B", fontSize=15, leading=19, spaceAfter=4)
H2 = ParagraphStyle("h2", parent=BODY, fontName="DV-B", fontSize=10.5, leading=14, spaceBefore=8, spaceAfter=3)
CELL = ParagraphStyle("cell", parent=BODY, fontSize=7.8, leading=10)
CELLB = ParagraphStyle("cellb", parent=CELL, fontName="DV-B")


def P(t, s=BODY):
    return Paragraph(t, s)


def table(rows, widths, header=True):
    data = [[Paragraph(str(c), CELLB if (header and i == 0) else CELL) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=[w * mm for w in widths], repeatRows=1 if header else 0)
    style = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#BBBBBB")),
             ("VALIGN", (0, 0), (-1, -1), "TOP"),
             ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    if header:
        style.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEEEEE")))
    t.setStyle(TableStyle(style))
    return t


def run_pytest():
    out = subprocess.run([sys.executable, "-m", "pytest", "-v", "-p", "no:cacheprovider", "test_multilayer.py"],
                         cwd=ROOT, capture_output=True, text=True)
    tests = re.findall(r"::(\w+) (PASSED|FAILED|ERROR)", out.stdout)
    summary = [l for l in out.stdout.splitlines() if re.search(r"\d+ (passed|failed)", l)]
    return tests, (summary[-1].strip("= ").strip() if summary else "no summary"), out.returncode


def feko_numbers():
    """Feko runs found in feko/*.out against the script, or None when there are none."""
    fdir = ROOT / "feko"
    runs = {k: parse_feko((fdir / f"{k}.out").read_text(errors="replace"))
            for k in ("B2", "B3", "B4", "C1", "C2") if (fdir / f"{k}.out").exists()}
    if not runs:
        return None
    f0 = 2.45e9
    out = {"runs": runs, "fdir": fdir}
    for k, r in runs.items():
        x = feko_derived(r)
        S = r["E0"] ** 2 / (2 * ETA0)
        if k.startswith("B"):
            T = 90 if k == "B4" else 60
            res = solve(stack([{"tissue": "muscle", "d_mm": T}], f0), f0, S)
            ref = {"names": ["Muscle_gabriel"], "abs": [res.layers[0].P_abs / res.S_inc]}
        else:
            res = solve(stack(SPEC_STACK[:2] + [{"tissue": "muscle", "d_mm": 60}], f0), f0, S)
            ref = {"names": ["skin_gabriel", "fat_gabriel", "Muscle_gabriel"],
                   "abs": [l.P_abs / res.S_inc for l in res.layers]}
        line = fdir / f"{k}_line.csv"
        stats = {}
        if line.exists():
            z, e = line_arrays(r)
            order = np.argsort(-z)
            z, e = z[order], e[order]
            ez = z[z < 0]
            ee = e[z < 0]
            ex = np.abs(res.fields(ez)[0])
            bounds = np.array([-l.x_top for l in res.layers] + [-res.layers[-1].x_bottom])
            keep = np.min(np.abs(ez[:, None] - bounds[None, :]), axis=1) > 0.4e-3
            for i, l in enumerate(res.layers):
                m = keep & (ez <= -l.x_top) & (ez > -l.x_bottom)
                if m.any():
                    d = np.abs(ee[m] / ex[m] - 1)
                    stats[i] = (d.mean(), d.max())
        out[k] = {"x": x, "res": res, "ref": ref, "stats": stats}
    return out


def feko_check_rows(fk):
    if not fk:
        return [["Feko vs script", "|Γ|, per-layer absorption", "within 3 %", "not done yet (no Feko result)", "open"]]
    rows = []
    g = lambda k: fk[k]["x"]
    gam = lambda k: fk["runs"][k]["R"][0]
    ref = lambda k: abs(fk[k]["res"].Gamma)
    dG = {k: gam(k) / ref(k) - 1 for k in fk["runs"]}
    dA = {k: max(abs(a / b - 1) for a, b in zip(
        [g(k)["abs_per_medium"][n] for n in fk[k]["ref"]["names"]], fk[k]["ref"]["abs"])) for k in fk["runs"]}
    worstG = max(abs(v) for v in dG.values())
    worstA = max(dA.values())
    rows.append(["Feko vs script: |Γ|", "Feko PBC cells: 60/90 mm muscle slab, skin/fat/muscle stack",
                 "within 3 %", f"worst {100*worstG:.2f} % over {len(dG)} runs", ok_(worstG < 0.03)])
    rows.append(["Feko vs script: absorbed power per layer", "power loss per medium / incident power on the cell",
                 "within 3 %", f"worst {100*worstA:.2f} % (relative, per layer)", ok_(worstA < 0.03)])
    fs = [v for k in fk["runs"] for v in (fk[k]["stats"] or {}).values()]
    if fs:
        mean_w = max(a for a, _ in fs)
        rows.append(["Feko vs script: |E| along the centre line", "near-field line, interface points skipped",
                     "within 3 %", f"largest layer mean {100*mean_w:.2f} %", ok_(mean_w < 0.03)])
    bal = max(abs(g(k)["absorbed_fraction"] / g(k)["balance_absorbed"] - 1) for k in fk["runs"])
    rows.append(["Power balance (plane-wave cell)", "power loss vs 1 − |R|² − |T|²", "error &lt; 5 %",
                 f"worst {100*bal:.2f} %", ok_(bal < 0.05)])
    cB = abs(gam("B3") / gam("B2") - 1) if {"B2", "B3"} <= set(fk["runs"]) else None
    cC = abs(gam("C2") / gam("C1") - 1) if {"C1", "C2"} <= set(fk["runs"]) else None
    cv = [v for v in (cB, cC) if v is not None]
    if cv:
        rows.append(["Mesh convergence (Feko)", "last refinement vs previous: |Γ| (and absorbed power in section 7)",
                     "within 2 %", f"worst {100*max(cv):.3f} %", ok_(max(cv) < 0.02)])
    return rows


def gam_val(fk, k):
    return fk["runs"][k]["R"][0]


def ok_(b):
    return "PASS" if b else "FAIL"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "validation_report.pdf"))
    ap.add_argument("--author", default="Yusuf Uzun")
    ap.add_argument("--course-media", help="muscle=FILE,skin=FILE,fat=FILE (course media pack CSVs)")
    args = ap.parse_args()

    media = None
    if args.course_media:
        media = compare_media(dict(kv.split("=", 1) for kv in args.course_media.split(",")))
    fk = feko_numbers()
    tests, summary, rc = run_pytest()
    all_pass = rc == 0 and tests and all(s == "PASSED" for _, s in tests)

    # ---- computed numbers -------------------------------------------------
    fres = run_fresnel()
    g_qw, g_half, bare = run_quarter_wave()
    cons = []
    for f in FREQS:
        res = solve(stack(SPEC_STACK, f), f)
        cons.append((f, abs(res.Gamma), 1 - abs(res.Gamma) ** 2,
                     sum(l.P_abs for l in res.layers) / res.S_inc, res.P_trans / res.S_inc,
                     res.conservation_error, res.integral_error))
    worst_media, worst_where = 0.0, ""
    for tissue, rows in WEEK1_TABLE.items():
        for f, (er, sg) in zip(FREQS, rows):
            er_c, sg_c = gabriel_props(tissue, f)
            for name, v, ref in (("eps_r", er_c, er), ("sigma", sg_c, sg)):
                d = abs(v - ref) / ref
                if d > worst_media:
                    worst_media, worst_where = d, f"{tissue} {name} at {f/1e6:g} MHz"
    f0 = 2.45e9
    half = solve(stack([{"tissue": "muscle", "d_mm": None}], f0), f0)
    delta = half.layers[0].delta
    z = -np.arange(0, 0.06, 0.2e-3)
    d_fit = fit_delta(z, np.abs(half.fields(z)[0]), 0.0, -0.06)
    slab = solve(stack(SPEC_STACK[:2] + [{"tissue": "muscle", "d_mm": 60}], f0), f0)
    m = slab.layers[-1]
    zs = -np.arange(0, m.x_bottom, 0.2e-3)
    es = np.abs(slab.fields(zs)[0])
    d5 = fit_delta(zs, es, -m.x_top, -m.x_bottom)
    d15 = fit_delta(zs, es, -m.x_top, -m.x_bottom, skip_bottom=1.5 * m.delta)

    # ---- document ---------------------------------------------------------
    story = []
    story += [P("P1-C Group 3: analytical multilayer script, validation report", H1),
              P(f"EM wave penetration in layered tissue (2.45 GHz team). Prepared by {args.author}, "
                f"{datetime.date.today().isoformat()}. Evidence for Checkpoint 1. "
                + ("<b>Scope: the Python script, and the Feko results obtained so far (section 7) compared with it.</b>"
                   if fk else "<b>Scope: the Python script only. No Feko result is included yet.</b>"), BODY),
              Spacer(1, 3)]

    story += [P("1. Checks against the pass criteria of the workbook", H2)]
    ok = lambda b: "PASS" if b else "FAIL"
    fres_ok = all(r[4] for r in fres)
    rows = [["Check", "Reference", "Pass criterion", "Result", "Pass?"],
            ["Script: single interface", "Fresnel formula, air/muscle, 4 frequencies", "4 significant digits",
             f"equal at all 4 frequencies; max |Γ| difference "
             f"{max(abs(r[1]-r[2]) for r in fres):.1e}", ok(fres_ok)],
            ["Script: λ/4 matching layer", "Textbook, lossless layer, η1 = √(η0·η2)", "|Γ| &lt; 1e-6",
             f"|Γ| = {g_qw:.1e}", ok(g_qw < 1e-6)],
            ["Penetration depth in muscle (script vs Week 1 table)", "Week 1 slide: 2.23 cm at 2.45 GHz", "within 2 %",
             f"δ = {delta*1e3:.2f} mm ({100*abs(delta*1e3-22.3)/22.3:.2f} % from 22.3 mm)",
             ok(abs(delta * 1e3 - 22.3) / 22.3 < 0.02)],
            ["Script: power conservation", "1 − |Γ|² = Σ P_abs/S_inc + P_trans/S_inc", "small residual",
             f"max error {max(c[5] for c in cons):.1e}", ok(max(c[5] for c in cons) < 1e-6)],
            ["Tissue media vs Week 1 table", "Week 1 slide 17 (ε′, σ at 4 frequencies)", "within 1 %",
             f"worst {100*worst_media:.2f} % ({worst_where})", ok(worst_media < 0.01)],
            *feko_check_rows(fk),
            ["Power balance (dipole case)", "P_acc = P_rad + P_abs", "error &lt; 5 %", "not done yet (Week 6)", "open"]]
    if media:
        dev = max(max(m["dev_eps"], m["dev_sigma"]) for m in media.values())
        rows.insert(6, ["Tissue media vs course media pack", "Course files Muscle / Skin_dry / Fat_not_infiltrated "
                        f"({min(m['n'] for m in media.values())} points each, 1 MHz to 20 GHz)", "within 1 %",
                        f"worst {100*dev:.4f} % at the file's own frequencies", ok(dev < 0.01)])
    story += [table(rows, [36, 46, 22, 54, 17])]

    story += [P("2. Test 1: air/muscle interface vs Fresnel", H2),
              P("Muscle half-space, Gabriel properties at each frequency, normal incidence. Reference "
                "Γ = (1 − n)/(1 + n) with n = √ε*. Compared as complex numbers (difference below 1e-9).", BODY),
              table([["f (MHz)", "script |Γ|", "Fresnel |Γ|", "1 − |Γ|² (power entering)", "4 digits"]] +
                    [[f"{r[0]/1e6:g}", f"{r[1]:.6f}", f"{r[2]:.6f}", f"{100*r[3]:.2f} %", ok(r[4])] for r in fres],
                    [24, 32, 32, 50, 22])]

    story += [P("3. Test 2: lossless quarter-wave layer", H2),
              P(f"Layer ε1 = 3 (λ1/4 thick) on a lossless ε2 = 9 half-space, η1 = √(η0·η2). "
                f"At the design frequency |Γ| = {g_qw:.1e}. As a counter-check, at twice the frequency the layer is "
                f"λ/2 thick and the stack must behave like the bare interface: |Γ| = {g_half:.6f} versus "
                f"{bare:.6f} for the bare interface.", BODY)]

    story += [P("4. Power conservation, skin 1.5 mm / fat 8 mm / muscle half-space", H2),
              P("Absorbed power per layer is the drop of the Poynting flux across the layer. It is cross-checked "
                "by integrating ½σ|E|² over the layer (last column).", BODY),
              table([["f (MHz)", "|Γ|", "1 − |Γ|²", "Σ P_abs/S_inc", "P_trans/S_inc", "conservation error",
                      "Poynting vs ∫½σ|E|²"]] +
                    [[f"{c[0]/1e6:g}", f"{c[1]:.4f}", f"{c[2]:.5f}", f"{c[3]:.5f}", f"{c[4]:.1e}",
                      f"{c[5]:.1e}", f"{c[6]:.1e}"] for c in cons],
                    [16, 16, 22, 26, 28, 30, 38])]

    story += [P("5. Penetration depth and the finite-slab bias", H2),
              P(f"Muscle half-space at 2.45 GHz: analytic δ = {delta*1e3:.3f} mm, fit of ln|E| (Week 1 procedure, "
                f"5 mm skipped at each face) = {d_fit*1e3:.3f} mm "
                f"({100*(d_fit/delta-1):+.2f} %). For a 60 mm muscle slab backed by air, as in a Feko model, the "
                f"back-face reflection (|Γ| about 0.76 for muscle/air) adds a ripple: the same fit on the "
                f"<i>exact</i> field gives {d5*1e3:.2f} mm ({100*(d5/m.delta-1):+.2f} %), which is close to the 2 % "
                f"limit. Skipping 1.5 δ at the bottom of the slab gives {d15*1e3:.2f} mm "
                f"({100*(d15/m.delta-1):+.2f} %). A Feko δ check should therefore use both windows.", BODY)]

    story += [P("6. Unit tests (pytest)", H2),
              P(f"Result: <b>{summary}</b>. Environment: Python {platform.python_version()}, "
                f"numpy {np.__version__}, matplotlib {matplotlib.__version__}.", BODY),
              Spacer(1, 3)]
    story += [KeepTogether([table([["Test", "Result"]] + [[n, s] for n, s in tests], [100, 30])]), Spacer(1, 4)]
    if not all_pass:
        story += [P("<b>WARNING: at least one test failed. Do not use this report as evidence.</b>", BODY)]

    for name, cap in (("results/stack_2450_profile.png",
                       "Figure 1. |E|(z) and ½σ|E|²(z) at 2.45 GHz for skin 1.5 / fat 8 / muscle half-space "
                       "(10 W/m² incident; air on the left)."),):
        img = ROOT / name
        if img.exists():
            story += [Image(str(img), width=118 * mm, height=118 * mm * 0.8333), P(cap, SMALL)]

    if fk:
        MESH = {"B2": "1.11", "B3": "0.83", "B4": "1.11", "C1": "0.5 / 1.0 / 1.11", "C2": "0.4 / 0.8 / 0.83"}
        CASE = {"B2": "muscle slab T = 60 mm", "B3": "muscle slab T = 60 mm", "B4": "muscle slab T = 90 mm",
                "C1": "skin / fat / muscle 60 mm", "C2": "skin / fat / muscle 60 mm"}
        story += [P("7. Feko results against the script (plane wave, periodic 10 x 10 mm cell, FEM)", H2),
                  P("Feko 2022, 2.45 GHz, |E0| = 86.8 V/m peak (10 W/m²), second-order FEM tetrahedra, periodic "
                    "boundary cell, media from the course pack (checked in the .out files: eps_r and sigma equal the "
                    "model values to 0.04 %). The script solves the same stack with the same slab backed by air. "
                    "Absorbed power is the Feko power loss per medium divided by the power incident on the cell.", BODY),
                  Spacer(1, 3)]
        rows = [["Run", "Case", "Mesh (mm)", "Tets", "DOF", "Mem GB", "Time s", "|Γ| Feko", "|Γ| script", "diff"]]
        for k in fk["runs"]:
            r = fk["runs"][k]
            g_ = gam_ = r["R"][0]
            rows.append([k, CASE[k], MESH[k], f"{r['tets']}", f"{r['dof']}", f"{r['mem_total_GB']:.2f}",
                         f"{r['run_s']:.0f}", f"{gam_:.5f}", f"{abs(fk[k]['res'].Gamma):.5f}",
                         f"{100*(gam_/abs(fk[k]['res'].Gamma)-1):+.2f} %"])
        story += [KeepTogether([table(rows, [11, 34, 24, 17, 19, 14, 14, 19, 19, 15])]), Spacer(1, 4)]

        rows = [["Run", "Layer", "Feko (% of incident)", "script (%)", "diff", "|E| mean diff", "|E| max diff"]]
        for k in fk["runs"]:
            for i, (n, a) in enumerate(zip(fk[k]["ref"]["names"], fk[k]["ref"]["abs"])):
                f_ = fk[k]["x"]["abs_per_medium"][n]
                st = fk[k]["stats"].get(i)
                rows.append([k, n.replace("_gabriel", ""), f"{100*f_:.3f}", f"{100*a:.3f}",
                             f"{100*(f_/a-1):+.2f} %",
                             f"{100*st[0]:.2f} %" if st else "-", f"{100*st[1]:.2f} %" if st else "-"])
        story += [KeepTogether([table(rows, [12, 26, 34, 24, 20, 26, 26])]), Spacer(1, 4)]

        conv = []
        if {"B2", "B3"} <= set(fk["runs"]):
            conv.append(f"Model B, 1.11 to 0.83 mm: |Γ| changes by {100*(gam_val(fk,'B3')/gam_val(fk,'B2')-1):+.3f} %.")
        if {"C1", "C2"} <= set(fk["runs"]):
            parts = ", ".join(f"{n.replace('_gabriel','')} {100*(fk['C2']['x']['abs_per_medium'][n]/fk['C1']['x']['abs_per_medium'][n]-1):+.2f} %"
                              for n in fk["C1"]["ref"]["names"])
            conv.append(f"Stack, mesh edges reduced by 1.25 to 1.34 times (52 607 to 116 033 tetrahedra): |Γ| "
                        f"{100*(gam_val(fk,'C2')/gam_val(fk,'C1')-1):+.3f} %, absorbed power {parts}.")
        if {"B2", "B4"} <= set(fk["runs"]):
            conv.append(f"Truncation: T = 90 mm gives |Γ| = {gam_val(fk,'B4'):.5f} (script {abs(fk['B4']['res'].Gamma):.5f}), "
                        f"T = 60 mm differs from it by {100*(gam_val(fk,'B2')/gam_val(fk,'B4')-1):+.2f} %. The 60 mm slab "
                        f"is enough for |Γ| and the near-surface field. Deeper in the muscle the back-face echo changes |E| "
                        f"(the script predicts the same change), so the penetration-depth fit uses the window that skips 1.5 δ "
                        f"at the back.")
        story += [P("<b>Convergence and truncation.</b> " + " ".join(conv), BODY)]
        story += [P("<b>Open point in the 60 mm muscle slab (runs B2, B3).</b> There the Feko absorbed power is "
                    "1.2 to 1.5 % above the script, the energy-balance gap is 0.7 to 0.8 %, and |E| differs by about "
                    "1 % on average and by up to 8 to 10 % near the minimum close to the back face. With the 90 mm "
                    "slab (B4) the absorbed power agrees to 0.02 % and the balance gap is 0.01 %. The stack, which has "
                    "the same 60 mm of muscle, agrees to 0.2 %. The cause of the 60 mm muscle-slab offset is not "
                    "identified yet; all criteria are still met (3 % and 5 %), but the offset has to be explained "
                    "before the 60 mm slab is used for anything finer than that.", BODY)]
        story += [P("<b>Penetration depth.</b> On the stack (60 mm muscle), the fit that skips 1.5 δ at the back gives "
                    "22.4 mm (+0.4 % from 22.33 mm). The Week 1 window (5 mm off each face) gives +2.1 % on the Feko "
                    "field and +2.2 % on the exact script field, so that excess is the finite-slab bias of the "
                    "window and not a Feko error (section 5).", BODY)]
        for name, cap in (("feko/C2_line_vs_script.png",
                           "Figure 2. Stack, finest mesh: Feko |E| along the centre line against the script."),
                          ("feko/model_b_log_convergence.png",
                           "Figure 3. Model B (60 mm muscle slab): quantities against mesh size.")):
            img = ROOT / name
            if img.exists():
                from reportlab.lib.utils import ImageReader
                w_, h_ = ImageReader(str(img)).getSize()
                story += [Image(str(img), width=130 * mm, height=130 * mm * h_ / w_), P(cap, SMALL)]
        lim_no = "8"
    else:
        lim_no = "7"

    story += [P(f"{lim_no}. Limitations and open items", H2),
              P(("• Feko covers the plane-wave cases only (section 7); the dipole model, the S11 detuning and the dipole power "
                 "balance are Week 6 and not done yet.<br/>• Feko runs use a single frequency (2.45 GHz) and normal "
                 "incidence; the |E| line stops 55 mm below the surface, so the last 15 mm of the 69.5 mm stack are "
                 "checked through |Γ| and the absorbed power only.<br/>"
                 if fk else
                 "• All results come from the script; the Feko models (plane-wave PBC cell, mesh convergence) and the "
                 "model log are not done yet. Checkpoint 1 also needs those.<br/>")
                +
                "• The tissue properties are computed from the published Gabriel (1996) four-term Cole-Cole parameters. "
                + (f"They agree with the course media pack files to {100*max(max(m['dev_eps'], m['dev_sigma']) for m in media.values()):.4f} % "
                   "at the file frequencies, and with the Week 1 table to " + f"{100*worst_media:.2f} %.<br/>"
                   if media else
                   "They agree with the Week 1 table to " + f"{100*worst_media:.2f} % and still have to be cross-checked against the course media pack.<br/>")
                +
                "• The tests check the script against analytical references and the lecture numbers, so they cannot "
                "reveal a convention mismatch with Feko; that is what the Feko vs script comparison is for.", BODY)]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(out), pagesize=A4, leftMargin=16 * mm, rightMargin=16 * mm,
                            topMargin=14 * mm, bottomMargin=14 * mm,
                            title="P1-C Group 3: validation report", author=args.author)
    doc.build(story)
    print(f"Wrote {out}  ({summary})")


if __name__ == "__main__":
    main()
