# P1-C prediction sheet — Group 3 (2.45 GHz)

Name: ______________________    Role this week: ______________    Date: __________

Fill this in **alone, on paper, before** running the script or Feko. Give a number, a sketch and a
one-line reason for each item. It is graded for completion, not correctness, and comes back in
Week 6 for your reflection. Do not open `prediction_answers.md` until the team has compared sheets.

Stack: air / skin 1.5 mm / fat 8 mm / muscle (half-space or 60 mm), plane wave at normal incidence,
10 W/m² incident (|E0| = 86.8 V/m peak).

## A. Single interface (Week 1 recap)

1. |Γ| for air → muscle at 2.45 GHz: ________   Reason: ______________________________
2. Fraction of incident power entering muscle: ________ %
3. Penetration depth δ in muscle (field 1/e): ________ mm. At what depth has the *power* dropped to 1/e? ________ mm
4. Wavelength in muscle λ_t: ________ mm → mesh edge for λ_t/10: ________ mm

## B. Layered stack (Sprint 1–2)

5. |Γ| for the full stack: larger / smaller / about the same as bare muscle?  ________   Number: ________
6. Rank the layers by absorbed power (most → least): ________ > ________ > ________
7. Percentage of the power that enters the tissue and is absorbed in **skin** (only 1.5 mm thick): ________ %
8. Sketch |E| vs z from z = +50 mm (air) to −60 mm. Mark the interfaces and where |E| jumps or changes slope.

```
|E|
 ^
 |
 |
 |
 +-----------------------------------------> z (mm)
 +50          0    -1.5   -9.5                  -60
     air      skin  fat        muscle
```

9. Sketch ½σ|E|² vs z on the same axis. In which layer is it largest per unit volume? ________
10. Is |E| continuous at each interface? Is ½σ|E|² continuous? Why? ______________________________

## C. Fat-thickness study (T5–T6 sweep, our explanation)

11. Fat thickness that gives the **most** power in muscle at 2.45 GHz: ________ mm
12. Is there a fat thickness where muscle absorbs **more** than with no fat at all? Yes / No. Range: ________ mm
13. Period of the pattern in fat thickness (hint: λ_fat): ________ mm

## D. Across frequency (class jigsaw)

14. Which frequency gives the smallest |Γ| for the stack: 403 / 915 / 2450 / 5800 MHz? ________
15. At 5.8 GHz, which layer absorbs the most? ________

## E. Dipole 10 mm above the stack (Sprint 3)

16. Compared with free space, the resonant frequency over tissue will go up / down / stay. By about ________ %
17. Radiation efficiency e_rad over tissue: ________ % (free-space dipole ≈ 100 %)
18. Near the surface, is more or less of the absorbed power in skin + fat than for the plane wave? ________  Why? ______________________________

**The number I am least sure about, and why:** ____________________________________________
