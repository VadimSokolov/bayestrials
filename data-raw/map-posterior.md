# Mapping: ch05_helpers.py, ch05_sequential.py -> R/posterior.R
- ch05_helpers.z_test_one_sided -> G1 `z_test_binary` (basics.R), not re-implemented here
- ch05_helpers.prob_superior -> core `prob_superior(a1, b1, a0, b0, method = "mc")` with `beta_post()` for the posterior parameters
- ch05_helpers.prob_superior_vec -> core `prob_superior` (vectorized over arguments; `method = "mc"` or "exact")
- ch05_helpers.hpd_beta -> `hpd_beta` (returns named c(lower, upper))
- ch05_helpers.normal_update -> core `normal_update` (returns list(mean, sd))
- ch05_sequential.looks_data, pr_superior, pooled_z -> internal to `sim_sequential_looks` (dropped as separate functions)
- ch05_sequential.run -> `sim_sequential_looks` (returns a list with the same names as the Python dict)
