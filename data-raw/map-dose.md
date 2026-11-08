# ch15_dose.py -> R/dose.R
- `TARGET`, `SCENARIOS`, `TOXIC` -> argument defaults (`target = 0.25`, `toxic = 0.35`) and `dose_scenarios()`
- `true_mtd(p, target)` -> `dose_true_mtd(p, target)`
- `three_plus_three(p, rng)` -> `dose_three_plus_three(p)` (global RNG)
- `crm(p, rng, skeleton, N, cohort)` -> `dose_crm(p, skeleton, N, cohort, target)` (grid and prior built inside)
- `GRID`, `PRIOR` -> dropped: built inside `dose_crm`
- `boin_bounds(target, f1, f2)` -> `dose_boin_bounds()` (named vector `escalate`, `deescalate`; tested against BOIN::get.boundary)
- `pava(r, w)` -> `dose_pava(r, w)`
- `boin(p, rng, N, cohort)` -> `dose_boin(p, N, cohort, target)`
- `run(design, p, reps, seed)` -> `dose_simulate(design, p, reps, seed)`; `design` is a function of `p`; `sel[1]` is "no dose", then doses 1..K
