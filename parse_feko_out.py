"""Read a Feko .out file (plane wave on a periodic cell) and report what Checkpoint 1 needs.

Usage:
  python parse_feko_out.py ModelB.out
  python parse_feko_out.py ModelB.out --run B1 --change "baseline, mesh 1.67 mm" --mesh 1.67 \
        --line-out feko/B1_line.csv --append feko/model_b_log.csv

Prints: medium values at the solution frequency, source amplitude, tetrahedra and FEM degrees of
freedom, reflection and transmission coefficient, power loss per medium as a share of the incident
power on the cell, energy balance, memory and times. With --line-out it writes the near-field line
(z in metres, total |E| peak) for model_log.py / compare_feko.py. With --append it adds one row to
the model log CSV (run, date, change, solver, mesh_mm, unknowns, mem_GB, time_min, gamma, line_file, notes).

Caveat: if the output contains "READING SOLUTION COEFFICIENTS FROM FILE" the solver re-used an earlier
solution, so its times and memory are not those of a fresh solve; the script says so.
"""
import argparse
import csv
import re
from pathlib import Path

import numpy as np

from tissue import ETA0

FLOAT = r"[-+]?\d+\.?\d*(?:[Ee][-+]?\d+)?"


def _section(text, start, ends):
    i = text.find(start)
    if i < 0:
        return ""
    j = min([k for k in (text.find(e, i + len(start)) for e in ends) if k >= 0] or [len(text)])
    return text[i:j]


def _to_gb(value, unit):
    return float(value) * {"kbyte": 1 / 1024 ** 2, "mbyte": 1 / 1024, "gbyte": 1.0}[unit.lower()]


def parse(text):
    d = {}
    g = lambda pat: re.search(pat, text)

    m = g(r"Number of degrees of freedom for FEM:\s*(\d+)")
    d["dof"] = int(m.group(1)) if m else None
    first = g(r"Number of first order tetrahedra:\s*(\d+)")
    second = g(r"Number of second order tetrahedra:\s*(\d+)")
    d["tets"] = (int(first.group(1)) if first else 0) + (int(second.group(1)) if second else 0)

    e0 = [float(g(rf"\|E0{c}\|\s*=\s*({FLOAT})").group(1)) if g(rf"\|E0{c}\|\s*=\s*({FLOAT})") else 0.0
          for c in "XYZ"]
    d["E0"] = float(np.sqrt(sum(x * x for x in e0)))
    vec = {}
    for name in ("A1X", "A1Y", "A1Z", "A2X", "A2Y", "A2Z"):
        m = g(rf"{name}\s*=\s*({FLOAT})")
        vec[name] = float(m.group(1)) if m else None
    if all(v is not None for v in vec.values()):
        a1 = np.array([vec["A1X"], vec["A1Y"], vec["A1Z"]])
        a2 = np.array([vec["A2X"], vec["A2Y"], vec["A2Z"]])
        d["cell_area"] = float(np.linalg.norm(np.cross(a1, a2)))
    else:
        d["cell_area"] = None
    f = g(r"Frequency in Hz:\s*FREQ\s*=\s*(" + FLOAT + ")")
    d["f_Hz"] = float(f.group(1)) if f else None

    media = []
    sec = _section(text, "DATA FOR DIELECTRIC MEDIA", ["EXCITATION", "DATA FOR THE GREEN"])
    for line in sec.splitlines():
        p = line.split()
        if len(p) >= 10 and p[0].isdigit():
            try:
                media.append({"index": int(p[0]), "eps_r": float(p[1]), "sigma": float(p[3]),
                              "tan_delta": float(p[4]), "wavelength_m": float(p[6]),
                              "density": float(p[7]), "name": " ".join(p[9:])})
            except ValueError:
                pass
    d["media"] = media

    losses = {}
    sec = _section(text, "POWER BUDGET IN THE DIELECTRIC", ["SUMMARY OF LOSSES"])
    for mm in media:
        mt = re.search(r"^\s*" + re.escape(mm["name"]) + r"\s+(.*)$", sec, re.M)
        if mt:
            last = mt.group(1).split()[-1]
            try:
                losses[mm["name"]] = float(last)
            except ValueError:
                pass
    d["loss_W"] = losses
    s = g(r"Sum of all losses:\s*(" + FLOAT + r")\s*W")
    d["loss_sum_W"] = float(s.group(1)) if s else None

    sec = _section(text, "TRANSMISSION AND REFLECTION COEFFICIENTS", ["read from buffer", "SUMMARY OF REQUIRED"])
    for line in sec.splitlines():
        nums = re.findall(FLOAT, line)
        if len(nums) == 11 and "E" in line:
            v = [float(x) for x in nums]
            d["T"] = (v[3], v[4])
            d["T_cross"] = (v[5], v[6])
            d["R"] = (v[7], v[8])
            d["R_cross"] = (v[9], v[10])
            break

    sec = _section(text, "VALUES OF THE ELECTRIC FIELD STRENGTH", ["VALUES OF THE MAGNETIC FIELD STRENGTH",
                                                                   "Near field request with name"])
    rows = []
    for line in sec.splitlines():
        p = line.split()
        if len(p) == 10 and p[0].lstrip("-").isdigit():
            try:
                v = [float(x) for x in p[1:]]
            except ValueError:
                continue
            rows.append((int(p[0]), v[0], v[1], v[2], v[3], v[5], v[7]))
    d["near"] = rows  # (medium, x, y, z, |Ex|, |Ey|, |Ez|)

    d["reused_solution"] = "READING SOLUTION COEFFICIENTS FROM FILE" in text
    t = g(r"total times:\s+(" + FLOAT + r")\s+(" + FLOAT + ")")
    d["cpu_s"], d["run_s"] = (float(t.group(1)), float(t.group(2))) if t else (None, None)
    mem = g(r"Sum of the peak memory of all processes:\s*(" + FLOAT + r")\s*(\w+)")
    d["mem_total_GB"] = _to_gb(mem.group(1), mem.group(2)) if mem else None
    mem1 = g(r"Peak memory usage during the whole solution:\s*(" + FLOAT + r")\s*(\w+)")
    d["mem_master_GB"] = _to_gb(mem1.group(1), mem1.group(2)) if mem1 else None
    its = re.findall(r"^\s*(\d+) of (\d+)\s+(" + FLOAT + r")\s*$", text, re.M)
    d["iterations"] = (int(its[-1][0]), float(its[-1][2])) if its else None
    return d


def derived(d):
    out = {}
    if d.get("E0"):
        out["S_inc"] = d["E0"] ** 2 / (2 * ETA0)
        if d.get("cell_area"):
            out["P_inc_W"] = out["S_inc"] * d["cell_area"]
            tot = sum(d["loss_W"].values())
            out["absorbed_fraction"] = tot / out["P_inc_W"]
            out["abs_per_medium"] = {k: v / out["P_inc_W"] for k, v in d["loss_W"].items()}
    if d.get("R"):
        r2 = d["R"][0] ** 2
        t2 = d["T"][0] ** 2 if d.get("T") else 0.0
        out["balance_absorbed"] = 1 - r2 - t2
        out["power_in"] = 1 - r2
    return out


def line_arrays(d):
    near = np.array([r for r in d["near"] if abs(r[1]) < 1e-12 and abs(r[2]) < 1e-12])
    z = near[:, 3]
    e = np.sqrt(near[:, 4] ** 2 + near[:, 5] ** 2 + near[:, 6] ** 2)
    return z, e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--run")
    ap.add_argument("--change", default="")
    ap.add_argument("--mesh", help="mesh edge in mm")
    ap.add_argument("--solver", default="Feko FEM (PBC)")
    ap.add_argument("--date", default="")
    ap.add_argument("--line-out")
    ap.add_argument("--append")
    args = ap.parse_args()

    d = parse(Path(args.out).read_text(errors="replace"))
    x = derived(d)

    print(f"Frequency {d['f_Hz']/1e9 if d['f_Hz'] else float('nan'):.4g} GHz, |E0| = {d['E0']:.4g} V/m, "
          f"cell area {d['cell_area']*1e6 if d['cell_area'] else float('nan'):.4g} mm^2")
    for m in d["media"]:
        print(f"  medium {m['index']} {m['name']}: eps_r = {m['eps_r']:.5f}, sigma = {m['sigma']:.5f} S/m, "
              f"tan delta = {m['tan_delta']:.5f}")
    print(f"Tetrahedra {d['tets']}, FEM degrees of freedom {d['dof']}")
    if d.get("R"):
        print(f"Reflection  |R| = {d['R'][0]:.6f} at {d['R'][1]:.2f} deg (cross-pol {d['R_cross'][0]:.2e})")
        print(f"Transmission |T| = {d['T'][0]:.6f} at {d['T'][1]:.2f} deg")
        print(f"1 - |R|^2 = {100*x['power_in']:.3f} %   1 - |R|^2 - |T|^2 = {100*x['balance_absorbed']:.3f} %")
    if "absorbed_fraction" in x:
        print(f"Power loss per medium as a share of the power incident on the cell "
              f"({x['P_inc_W']:.4g} W):")
        for k, v in x["abs_per_medium"].items():
            print(f"  {k}: {100*v:.3f} % ({d['loss_W'][k]:.5g} W)")
        bal = abs(x["absorbed_fraction"] - x["balance_absorbed"]) / x["balance_absorbed"] if d.get("R") else None
        if bal is not None:
            print(f"  loss vs 1 - |R|^2 - |T|^2: {100*bal:.2f} % apart")
    if d["iterations"]:
        print(f"Iterative solver: {d['iterations'][0]} iterations, last residual {d['iterations'][1]:.3e}")
    print(f"Memory: peak {d['mem_master_GB']:.4g} GB (master), {d['mem_total_GB']:.4g} GB (all processes); "
          f"times: cpu {d['cpu_s']} s, run {d['run_s']} s")
    if d["reused_solution"]:
        print("NOTE: this output re-used an earlier solution (READING SOLUTION COEFFICIENTS FROM FILE): "
              "times and memory are not those of a fresh solve.")

    if args.line_out:
        z, e = line_arrays(d)
        order = np.argsort(-z)
        Path(args.line_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.line_out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["z_m", "abs_E_Vpm"])
            for zi, ei in zip(z[order], e[order]):
                w.writerow([f"{zi:.6e}", f"{ei:.6f}"])
        print(f"Wrote {len(z)} points to {args.line_out}")

    if args.append:
        path = Path(args.append)
        new = not path.exists()
        with open(path, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["run", "date", "change", "solver", "mesh_mm", "unknowns", "mem_GB", "time_min",
                            "gamma", "line_file", "notes"])
            note = "solution re-used from file: time/memory not representative" if d["reused_solution"] else ""
            w.writerow([args.run or "", args.date, args.change, args.solver, args.mesh or "",
                        d["dof"] or d["tets"], f"{d['mem_total_GB']:.4g}" if d["mem_total_GB"] else "",
                        f"{d['run_s']/60:.3g}" if d["run_s"] is not None else "",
                        f"{d['R'][0]:.6f}" if d.get("R") else "",
                        Path(args.line_out).name if args.line_out else "", note])
        print(f"Appended a row to {path}")


if __name__ == "__main__":
    main()
