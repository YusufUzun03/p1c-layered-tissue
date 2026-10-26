# Feko guide: Model B (Sprint 1) and the layered stack (Sprint 2), 2.45 GHz

Based on Week 1 Demos 1–3 and the Week 4 Brief C slides. Reference numbers below come from this repo's
script. **CADFEKO menu names are not from the slides**, so check them against your Feko version; the
settings (what to enter, which values) are from the slides and the script.

## 0. Targets (from `python model_log.py ...` reference, 10 W/m² or any E0)

| Quantity | Muscle slab T = 60 mm | T = 90 mm (truncation check) |
| --- | --- | --- |
| \|Γ\| | 0.76165 | 0.76246 (half-space, Fresnel) |
| δ (muscle, field 1/e) | 22.33 mm | 22.33 mm |
| Power through the back face | 0.0083 W/m² (0.08 %) | 0.0006 W/m² |
| \|E\| at 10 mm depth, E0 = 1 V/m | 0.1544 V/m | 0.1528 V/m |
| λ_t in muscle | 16.73 mm | |

T = 60 mm vs 90 mm changes |Γ| by 0.1 % in the script, so the slab thickness is not a risk. Pass criteria
(Week 1 and 4): |Γ| and δ within 2 % of the script, convergence < 2 % on fields, memory ≤ 12 GB, ≤ 15 min.

## 1. Variables (use them everywhere)

`f0 = 2.45e9`, `cell = 10` mm, `T = 60` mm, `z_air = 50` mm, `mesh = 1.67` mm, `E0` (1 V/m or 86.8 V/m;
write it in the log, `model_log.py --E0 1` takes it).

## 2. Media (Week 1 Demo 1)

1. `muscle_gabriel`: dielectric, definition **Frequency list (linear interpolation)**, permittivity from
   the course `Muscle_eps.csv`, conductivity from `Muscle_sigma.csv` (scale 1, comma delimiter), density 1090 kg/m³.
   Backup files from `python gabriel_media.py` are in `media/` (same format guess: frequency [Hz], value).
2. Enter **σ**, not tan δ. If you must use tan δ, write it in the log.
3. Optional check medium `muscle_const`: εr = 52.7, σ = 1.74 S/m. Both muscles must agree at 2.45 GHz.
4. After the first solve, read εr and σ at 2450 MHz in the `.out` file and compare with
   `media/course_frequencies_check.csv` (within 1 %).
5. Make `skin` (dry) and `fat` (not infiltrated) the same way now; you need them in Sprint 2.

## 3. Model B geometry and source (Week 1 Demo 2)

1. Unit cell 10 × 10 mm in x–y with a **periodic boundary condition**, solver **FEM** (MoM-SEP is the alternative).
2. Muscle slab from z = 0 down to z = −T, air above.
3. Plane wave, normal incidence (θ = 180°, travelling toward −z), amplitude E0.
4. Near-field line through the cell centre from z = +50 mm to z = −(T − 5) mm. The brief asks for a **0.2 mm**
   step (Week 1 used 0.5 mm). Do not place a sample exactly on the interface.
5. Requests: the near-field line, the reflection coefficient of the periodic structure, and power loss per
   medium (or integrate ½σ|E|² per layer afterwards). In the Feko example for periodic structures the
   transmission/reflection coefficients are requested for the incident plane wave; I could not confirm
   the exact request names for the power loss per medium, so look for it in the request list.

## 4. Mesh convergence (Week 4)

Pick one quantity Q first (use |Γ|; also |E| at a point 10 mm into the muscle and δ). Three meshes in the muscle:

| Run | Edge | Rough tet count* |
| --- | --- | --- |
| B1 | λ_t/10 = 1.67 mm | ≈ 11 000 |
| B2 | λ_t/15 = 1.11 mm | ≈ 37 000 |
| B3 | λ_t/20 = 0.83 mm | ≈ 89 000 |

*Order of magnitude only (cell volume / (0.118 · edge³)); Feko's actual number, with air and grading, is what
you log. Refine locally at the tissue surface before refining globally. Log unknowns, peak memory and time
for every run in `feko/model_b_log_template.csv` (copy it to `feko/model_b_log.csv`).

## 5. After each solve

1. POSTFEKO: plot |E| vs z, export the line to CSV (z in m, |E| peak). Check that |E| decays smoothly with no ripple
   from the back face.
2. Put the file name in the log's `line_file` column, |Γ| in `gamma`.
3. Run `python model_log.py feko/model_b_log.csv` (add `--E0 1` if you used 1 V/m).
   It prints convergence, agreement with the script and the laptop budget, and writes the workbook rows.
4. Truncation test: solve once more with T = 90 mm (log as B4, run `--T-mm 90`); results change < 2 %.

## 6. Extending to the stack (Sprint 2)

Skin 1.5 mm, fat 8 mm, muscle 60 mm. Script reference (`inputs/stack_2450_feko_slab.json`):

| Quantity | Value (10 W/m² incident) |
| --- | --- |
| \|Γ\| | 0.52229 |
| Absorbed in skin / fat / muscle | 3.167 / 0.5515 / 3.5466 W/m² |
| Through the back face | 0.0070 W/m² |

Mesh: at least 3 elements across the skin (edge ≤ 0.5 mm there), fat ≥ 2 elements, muscle as in section 4.
Feko power loss per medium is in W per cell: divide by the cell area 1e-4 m² to get W/m² and compare with the
table (within 3 %). Use `python compare_feko.py inputs/stack_2450_feko_slab.json line.csv --gamma ... --pabs ...`.

## 7. If the numbers do not match

| Symptom | Check |
| --- | --- |
| \|Γ\| off by a few % | σ or εr entered wrong at 2.45 GHz (read the `.out` file), tan δ entered as σ |
| |E| off by a constant factor | E0 mismatch (1 V/m vs 86.8 V/m), peak vs RMS |
| δ off by about 2 % only in the slab | back-face echo; the 5 mm-skip fit is biased, use the 1.5δ-skip fit that `compare_feko.py` and `model_log.py` print |
| Jumps in |E| at the interfaces | a sample sits exactly on the interface, or the line is not through the cell centre |
| Strong ripple in |E| deep in the muscle | slab too thin or the mesh too coarse at the bottom |
| Run over 15 min or 12 GB | coarser mesh far from the surface first; PBC cell is already small |
