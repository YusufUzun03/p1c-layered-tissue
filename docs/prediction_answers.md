# P1-C prediction sheet — reference answers (open after the team comparison)

Computed with this repo's script (Gabriel 1996 media, exact transmission-line solution) and the
Week 1–3 slides. Part E is qualitative until the Feko dipole runs are done.

| # | Answer | Source |
| --- | --- | --- |
| 1 | \|Γ\| = 0.762 (0.7625; Γ = −0.7621 + j0.0250) | `python test_multilayer.py`; Week 1 slide 21 |
| 2 | 41.9 % (1 − \|Γ\|²) | same |
| 3 | δ = 22.3 mm; power 1/e depth = δ/2 = 11.2 mm | Week 1 slide 10 |
| 4 | λ_t = 16.7 mm → λ_t/10 = 1.7 mm (λ_t/15 = 1.1, λ_t/20 = 0.8) | Week 1 slide 21 |
| 5 | Smaller: 0.523. The skin + fat stack acts partly as a matching section | `run.py inputs/stack_2450.json` |
| 6 | muscle 48.7 % > skin 43.7 % > fat 7.6 % (of the power entering) | same |
| 7 | 43.7 %: skin is lossy (σ = 1.46 S/m) and sits where \|E\| is largest | same |
| 8 | Standing wave in air (SWR ≈ 3.2); \|E\| is continuous at every interface, rises slightly through the fat (low η), then decays as e^(−x/δ) in muscle | `results/stack_2450_profile.png` |
| 9 | Largest per volume in skin (≈ 2100 W/m³ vs ≈ 300 W/m³ at the top of muscle) | same |
| 10 | Tangential E and H are continuous; ½σ\|E\|² jumps because σ jumps (fat σ = 0.105 vs muscle 1.74 S/m) | theory |
| 11 | ≈ 4 mm (5.43 W/m² in muscle) | `python sweeps.py` |
| 12 | Yes: about 0.5–7 mm gives more than the no-fat value of 3.84 W/m². The 8 mm spec stack gives 3.54, slightly less | same |
| 13 | λ_fat/2 ≈ 26.5 mm; \|Γ\| minimum near 5.5 mm, maximum near 20 mm | same |
| 14 | 2.45 GHz of the four (0.523); the continuous curve has its minimum near 1.85 GHz | `results/freq_sweep.png` |
| 15 | Skin 45.7 %, fat 26.1 %, muscle 28.2 %: skin and fat together take over 70 % | same |
| 16 | Down: high-ε tissue in the reactive near field (λ/2π ≈ 19 mm) loads the dipole. The size of the shift is a Feko result (Week 3 Demo 3 sweep) | Week 3 slides 18, 24 |
| 17 | Drops well below 100 %; Week 3 quotes 10–50 % on-body. Feko number to come | Week 3 slide 18 |
| 18 | More in skin + fat: near fields decay faster than e^(−x/δ) and have normal components, so absorption is pulled toward the surface | Week 2 slide 7, Week 3 |
