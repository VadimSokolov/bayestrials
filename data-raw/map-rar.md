# ch14_rar.py -> R/rar.R
- `alloc_probs(rule, S, F, n_done, N, rng, K, ...)` -> `rar_alloc_probs(rule, S, F, n_done, N, K, ...)` (no rng argument; uses global RNG)
- `run(p, N, rule, drift, reps, batch, seed, n_strata, fixed_ctrl)` -> `rar_simulate(...)` (returns a named list; `seed = NULL` default)
- `RULES2`, `RULES4` -> internal `rar_rules(K)`
- `two_arm(reps, N, p)` -> `rar_two_arm()` (data.frame, one row per rule)
- `drift_sweep(drifts, reps, N)` -> `rar_drift_sweep()` (data.frame, rows per rule and drift)
- `four_arm(reps, N, p, drift)` -> `rar_four_arm()` (list with data.frames `alt`, `null`)
- `md_table(header, rows)` -> dropped: formatting helper, use knitr::kable
