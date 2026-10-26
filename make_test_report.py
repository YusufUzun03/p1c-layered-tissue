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
from tissue import gabriel_props

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "validation_report.pdf"))
    ap.add_argument("--author", default="Yusuf Uzun")
    ap.add_argument("--course-media", help="muscle=FILE,skin=FILE,fat=FILE (course media pack CSVs)")
    args = ap.parse_args()

    media = None
    if args.course_media:
        media = compare_media(dict(kv.split("=", 1) for kv in args.course_media.split(",")))
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
                f"{datetime.date.today().isoformat()}. Evidence for Checkpoint 1 (script part). "
                f"<b>Scope: the Python script only. No Feko result is included yet.</b>", BODY),
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
            ["Feko vs script", "|Γ|, per-layer absorption", "within 3 %", "not done yet (no Feko result)", "open"],
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

    story += [P("7. Limitations and open items", H2),
              P("• All results come from the script; the Feko models (plane-wave PBC cell, mesh convergence) and the "
                "model log are not done yet. Checkpoint 1 also needs those.<br/>"
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
