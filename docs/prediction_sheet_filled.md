# P1-C prediction sheet, filled in with the reasoning (Group 3, 2.45 GHz)

Study copy. Every answer has the number and how it is obtained, so you can redo it by hand and write
your own sheet. Marks: **[hand]** = you can compute it with a calculator and the Week 1 table;
**[script]** = needs the full transmission-line calculation, taken from this repo's script;
**[guess]** = a reasoned estimate, nothing computed yet (Feko will settle it).

Setup: air / skin 1.5 mm / fat 8 mm / muscle, plane wave, normal incidence, 10 W/m², |E0| = 86.8 V/m peak
(E0 = √(2·η0·S) = √(2·377·10)).

## A. Single interface

**1. |Γ| air→muscle = 0.762** [hand]
Week 1 table: muscle at 2.45 GHz has εr = 52.7 and σ = 1.74 S/m.
- σ/(ωε0) = 1.74 / (2π·2.45e9·8.854e-12) = 1.74 / 0.1363 = 12.76, so ε* = 52.7 − j12.76 (tan δ = 0.24, as in the slides).
- n = √ε*: |ε*| = 54.25, angle −13.6°, so n = 7.365 ∠ −6.8° = 7.31 − j0.87.
- Γ = (1 − n)/(1 + n) = (−6.31 + j0.87)/(8.31 − j0.87), \|Γ\| = 6.37 / 8.36 = 0.762.

**2. Power entering muscle = 41.9 %** [hand]: 1 − \|Γ\|² = 1 − 0.581 = 0.419.

**3. δ = 22.3 mm; power falls to 1/e at δ/2 = 11.2 mm** [hand]
α = (ω/c)·\|Im n\| = (2π·2.45e9 / 3e8)·0.872 = 51.3·0.872 = 44.8 Np/m, so δ = 1/α = 22.3 mm. The field falls as e^(−x/δ),
the power (∝ \|E\|²) as e^(−2x/δ), so its 1/e depth is half of that.

**4. λ_t = 16.7 mm, mesh edge λ_t/10 = 1.67 mm** [hand]: λ0 = c/f = 122.4 mm, λ_t = λ0 / Re(n) = 122.4 / 7.31 = 16.7 mm.

## B. Layered stack

**5. \|Γ\| of the stack: smaller than bare muscle, about 0.52** [script]
Reasoning by hand is qualitative only: skin (εr 38) is close to muscle, fat (εr 5.3) is very different. Together they are a
thin two-layer section between air and muscle, which is a partial matching layer, so reflection goes down. The number
needs the impedance recursion (Z_in), which the script does: 0.523.

**6. Absorption ranking: muscle > skin > fat (49 %, 44 %, 8 % of the power entering)** [script for the numbers]
Fat has σ = 0.105 S/m, about 17 times less than muscle (1.74), so it absorbs little even though it is 8 mm thick. Skin is only
1.5 mm but has σ = 1.46 S/m and sits at the surface where the field is strong. Muscle wins because it is thick (δ = 22 mm).

**7. Skin takes about 44 % of what enters** [hand, once Γ is known]
- Surface field: \|E(0)\| = E0·\|1 + Γ\| = 86.8·0.607 = 52.7 V/m (with Γ = −0.452 − j0.263 from the script).
- Skin is thin, so \|E\| barely changes across it: P_skin ≈ ½·σ·\|E\|²·d = 0.5·1.464·52.7²·1.5e-3 = 3.05 W/m².
- Entering power = 10·(1 − 0.523²) = 7.26 W/m², so skin ≈ 42 % (script: 43.7 %).

**8. Sketch of \|E\|(z)**
- Air: standing wave with SWR = (1 + 0.52)/(1 − 0.52) ≈ 3.2. Maximum E0(1 + \|Γ\|) = 132 V/m, minimum E0(1 − \|Γ\|) = 41 V/m, repeating every λ0/2 = 61 mm.
- At z = 0 the field is 53 V/m. It stays continuous at every interface and shows a small bump across the fat.
- At the top of the muscle it is about 19 V/m, then it decays as e^(−x/22.3 mm).

**9. Sketch of ½σ\|E\|²(z): largest in skin** [hand/script]
Skin: ½·1.464·53² ≈ 2000 W/m³ (script peak 2150). Fat: tiny, σ is small and E is similar. Muscle top: ½·1.74·19² ≈ 310 W/m³, then it decays as e^(−2x/δ).

**10. \|E\| is continuous; ½σ\|E\|² is not** [theory]
At a planar interface the tangential E and H are continuous (Maxwell's boundary conditions), so \|E\| is continuous for normal incidence.
The loss density has σ as a factor, and σ jumps (fat 0.105 vs muscle 1.74 S/m), so it jumps with it.

## C. Fat-thickness study

**11. Best fat thickness for muscle absorption: about 4 mm** [script]; a by-hand estimate would say "a few mm, well below λ_fat/4 = 13 mm".
**12. Yes: about 0.5–7 mm gives more than no fat** [script]. Physical reason: a thin low-ε layer between skin and muscle lowers the reflection.
**13. Period ≈ 26.6 mm** [hand]: λ_fat = λ0/√εr,fat = 122.4/√5.28 = 53.3 mm, and the reflection pattern repeats every λ_fat/2 = 26.6 mm.
(Minimum of \|Γ\| at 5.5 mm, maximum at 20 mm in the script.)

## D. Across frequency

**14. Smallest \|Γ\| of the four frequencies: 2450 MHz (0.52)** [script]. The other values are 0.77 (403), 0.65 (915) and 0.87 (5800).
**15. At 5.8 GHz skin absorbs the most (46 %), then muscle (28 %) and fat (26 %)** [script]. Here δ in skin is only 8.6 mm and in muscle 7.5 mm, so
the thin skin is a large part of the penetration depth and the field reaches the muscle weaker. This agrees with the instructor's tip that absorption concentrates in skin above ~6 GHz.

## E. Dipole 10 mm above the stack (all [guess], to be replaced by Feko)

**16. Resonance goes down, by roughly 3–5 %** [guess]
The dipole is 10 mm from the surface, inside the reactive near field (λ0/2π = 19.5 mm). Tissue with εr 38–53 loads it like added capacitance, so f_res moves down (Week 3: "f_r decreases as g decreases"). The size of the shift is my guess.

**17. Radiation efficiency roughly 30 %** [guess]: Week 3 quotes 10–50 % for on-body antennas. At g = 10 mm this is in the middle of that range.

**18. More of the absorbed power in skin + fat than for the plane wave** [guess from theory]
Near fields decay faster than a plane wave and have a field component normal to the surface, so the power is deposited closer to the surface.

**Number I am least sure about:** the dipole efficiency in item 17: nothing in the lecture notes pins it down for a 10 mm gap, and it decides how much power the tissue takes.
