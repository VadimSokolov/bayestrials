# ch04_sim.py to R (randomization.R)

- `confounded_trial(n, rng, effect, assign)`: `sim_confounded_trial(n, effect, assign = "doctor"|"coin", seed)`, returns named vector `naive`, `health_gap`.
- `simple_assign`: `assign_simple(n, seed)`.
- `block_assign(n, rng, sizes)`: `assign_block(n, sizes, seed)`.
- `stratified_block_assign(strata, rng, sizes)`: `assign_stratified(strata, sizes, seed)`.
- `minimization_assign(X, rng, p_follow)`: `assign_minimization(X, p_follow, seed)`.
- `covariate_imbalance(arm, x)`: same name.
- `perm_pvalue(y, arm, rng, B, restricted_blocks)`: `perm_test_arm(y, arm, nsim, restricted_blocks, seed)`, returns named vector `obs`, `p`. Restricted blocks are 1-based index vectors in R.
- `arnold_s(points, L)`: `arnold_stochasticity(points, L)`, returns named vector `S`, `expected`.
- `compress_ratio(bits)`: same name. Uses `memCompress(type = "gzip")` (zlib, default level; Python uses level 9), so values differ slightly.
- `adjusted_vs_naive(n, rng, effect, bands)`: `sim_adjusted_vs_naive(n, effect, bands, seed)`, returns named vector `adjusted`, `naive`.
- RNG streams differ from numpy, so seeded draws are not identical; agree within Monte Carlo error.
