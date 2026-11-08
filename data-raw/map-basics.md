# ch03_sim.py to R (basics.R)

- `P0, P1, N, MU_DIFF, SD`: dropped, chapter constants; define in the chunk.
- `risk_diff(y1, y0)`: dropped, use `mean(y1) - mean(y0)`.
- `sim_binary_trials(rng, p0, p1, n, reps)`: `sim_two_arm_binary(p0, p1, n, nsim, seed)`, returns data.frame `x0`, `x1`, `diff`.
- `z_binary(x0, x1, n)`: `z_test_binary(x0, x1, n)`.
- `perm_test(rng, y1, y0, reps, stat)`: `perm_test(y1, y0, nsim, seed)`, returns list `obs`, `sims`, `p_one`, `p_two`. Unused `stat` argument dropped.
- `wald_ci(x0, x1, n, level)`: `wald_ci(x0, x1, n, level)`, returns data.frame `lower`, `upper`.
