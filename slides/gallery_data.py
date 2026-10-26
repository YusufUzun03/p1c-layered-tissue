"""Collect the numbers for the Week 4 gallery slide into gallery_data.json (run from the repo root)."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from multilayer import Layer, fit_delta, make_layer, solve  # noqa: E402
from tissue import complex_eps, gabriel_eps, wave_params  # noqa: E402

FREQS = [403e6, 915e6, 2450e6, 5800e6]
STACK = [{"tissue": "skin", "d_mm": 1.5}, {"tissue": "fat", "d_mm": 8}, {"tissue": "muscle", "d_mm": None}]

rows = []
for f in FREQS:
    bare = solve([make_layer({"tissue": "muscle", "d_mm": None}, f)], f)
    n = np.sqrt(gabriel_eps("muscle", f))
    stack = solve([make_layer(s, f) for s in STACK], f)
    rows.append({"f_MHz": f / 1e6, "fresnel": abs((1 - n) / (1 + n)), "script": abs(bare.Gamma),
                 "stack": abs(stack.Gamma)})

F = 2450e6
g = rows[2]
lam1 = wave_params(complex_eps(3.0, 0, F), F)[3]
qw = solve([Layer("m", lam1 / 4, complex_eps(3.0, 0, F), F), Layer("l", None, complex_eps(9.0, 0, F), F)], F)
cons = max(solve([make_layer(s, f) for s in STACK], f).conservation_error for f in FREQS)

slab = solve([make_layer({"tissue": "muscle", "d_mm": 60}, F)], F)
m = slab.layers[0]
z = -np.arange(0, m.x_bottom, 0.2e-3)
absE = np.abs(slab.fields(z)[0])
d5 = fit_delta(z, absE, 0, -m.x_bottom)
dsafe = fit_delta(z, absE, 0, -m.x_bottom, skip_bottom=1.5 * m.delta)

data = {
    "rows": rows,
    "gamma_2450": g["script"], "gamma_err": abs(g["script"] - g["fresnel"]) / g["fresnel"],
    "qw_gamma": abs(qw.Gamma), "conservation": cons,
    "delta_mm": m.delta * 1e3, "fit5_err": (d5 - m.delta) / m.delta, "fitsafe_err": (dsafe - m.delta) / m.delta,
}
out = Path(__file__).with_name("gallery_data.json")
out.write_text(json.dumps(data, indent=2))
print(json.dumps(data, indent=2))
