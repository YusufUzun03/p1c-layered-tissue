# P1-C multilayer script (Group 3, 2.45 GHz)

Exact transmission-line solution for a plane wave at normal incidence on N planar tissue layers.
Same geometry and conventions as Feko and the course: air for z > 0, wave from +z, tissue for z < 0,
e^{jωt}, peak phasors, 10 W/m² incident (|E0| = 86.8 V/m).

```
pip install -r requirements.txt
python test_multilayer.py                 # Checkpoint 1 validation report
pytest -q                                 # the same checks as unit tests
python run.py inputs/stack_2450.json      # one case: Gamma, |E|(z), absorbed power per layer
python sweeps.py                          # Gamma vs f, per-layer split at 4 frequencies, fat sweep
python model_log.py feko/model_b_log.csv  # Checkpoint 1: Model B mesh convergence (add --E0 1 if Feko used 1 V/m)
python dipole_reference.py                # Week 6: free-space dipole reference
python compare_feko.py inputs/stack_2450_feko_slab.json feko_line.csv --gamma 0.522 \
       --pabs skin=3.17,fat=0.55,muscle=3.55        # Week 5: Feko vs script
```

## Files

| File | What it does |
| --- | --- |
| `tissue.py` | Gabriel 1996 four-term Cole-Cole model (skin, fat, muscle, bone), k, η, δ, λ |
| `multilayer.py` | Solver: input-impedance recursion, fields E(z), H(z), Poynting flux, ½σ\|E\|², SAR, δ fit |
| `run.py` | Input file → summary JSON, profile CSV (class schema) and plot |
| `sweeps.py` | Outputs 1, 3 and 4 of the brief |
| `compare_feko.py` | Overlays a POSTFEKO near-field export and runs the Week 5 checks |
| `dipole_reference.py` | Free-space dipole Z_in(L/λ), resonance length (induced EMF): the Week 2 / Week 6 reference |
| `compare_dipole.py` | Week 6: power balance, e_rad, per-layer shares and on-axis profile, dipole vs plane wave |
| `gabriel_media.py` | Gabriel εr(f), σ(f) lists for Feko (backup to the course media pack) |
| `docs/feko_model_b_guide.md` | Step-by-step Model B and stack set-up with the reference numbers |
| `model_log.py` | Feko runs (log CSV + optional line exports) → convergence checks, script comparison, rows for the workbook's Model log |
| `feko/model_b_log_template.csv` | Copy and fill one row per Feko solve |
| `docs/prediction_sheet.md` | Individual prediction sheet (answers in `docs/prediction_answers.md`, open after comparing) |
| `slides/` | Week 4 gallery slide (`gallery_slide_T3.pptx` is included). Rebuilding it needs the pptx theme helper `apply_theme.js`, which is not in this repo: `python slides/gallery_data.py`, then `node slides/gallery_slide.js` |
| `test_multilayer.py` | Fresnel (4 digits), λ/4 match (\|Γ\| < 1e-6), power conservation, media vs Week 1 table, δ, continuity |
| `inputs/*.json` | Cases: spec stack with muscle half-space, 60 mm slab as in Feko, explicit εr/σ |

## Input file

```json
{
  "team": "T3", "case": "plane_wave", "f_MHz": 2450, "S_inc_Wpm2": 10,
  "dz_mm": 0.2, "air_mm": 50, "halfspace_depth_mm": 60,
  "layers": [
    {"name": "skin",   "tissue": "skin",   "d_mm": 1.5},
    {"name": "fat",    "tissue": "fat",    "d_mm": 8},
    {"name": "muscle", "tissue": "muscle", "d_mm": null}
  ]
}
```

Each layer takes either `tissue` (Gabriel model at `f_MHz`) or explicit `eps_r` and `sigma` (S/m).
`d_mm: null` is a half-space; a finite last layer is backed by air, like a Feko slab.

## Outputs

- `*_summary.json`: Γ (complex and magnitude), power entering, absorbed power per layer in W/m² and as a fraction of incident and entering power, the conservation errors, and warnings.
- `*_profile.csv`: `team, f_MHz, case, z_mm, layer, absE_Vpm, pabs_Wpm3`, every 0.2 mm (the class CSV schema).
- `*_profile.png`: |E|(z) and ½σ|E|²(z).

## How the checks work

- Absorbed power per layer is the drop in Poynting flux across the layer. It is cross-checked against the integral of ½σ|E|² (agreement ~1e-7).
- Power conservation: 1 − |Γ|² = Σ P_abs / S_inc + P_trans / S_inc (error ~1e-9).
- δ fit: a 60 mm muscle slab backed by air reflects strongly at its back face. Fitting ln|E| with only 5 mm skipped at each face is then off by about ±2 %, even on the exact field. `fit_delta(..., skip_bottom=1.5*delta)` brings this under 1 %. `compare_feko.py` prints both fits.

## Status (Group 3, P1-C, 2.45 GHz)

- Analytical script: done and validated (`python test_multilayer.py`, `pytest -q`, report in `docs/validation_report.pdf`).
- Feko: CADFEKO v2022 project started (muscle frequency-list medium imported). No Feko results yet; the Feko vs script comparison, the mesh-convergence study and the model log are open.
- `docs/feko_model_b_guide.md` lists the steps and the reference numbers to compare against.

## Notes

- The script, tools and documents were developed with the help of an AI assistant (Claude).
- Course material (lecture slides, the course media pack, the workbook) is not included. The tissue properties are computed from the published Gabriel (1996) four-term Cole-Cole parameters; `media/` holds them as frequency lists for Feko (the course pack is the primary source).
- `docs/prediction_sheet*.md` and `docs/prediction_answers.md` are study copies for the individual prediction sheet.
