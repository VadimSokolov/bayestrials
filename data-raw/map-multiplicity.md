# ch09_multiplicity.py to R

- `bonferroni` -> `mult_bonferroni`
- `holm` -> `mult_holm`
- `hochberg` -> `mult_hochberg`
- `graph_test` -> `mult_graph`
- `ULCER` -> dataset `ulcer` (data.frame: operation, n, recurrences; counts illustrative)
- `ulcer_logodds` -> `shrink_logodds(ulcer)` returns list(y, v)
- `eb_fit` -> `shrink_eb_fit(y, v)` returns list(mu, tau, shrunk, sd)
- `norm_post_prob` -> `shrink_post_prob(est, se, prior_sd)` returns list(prob, precision)
- `simulate_procedures` -> `mult_simulate` (data.frame output; arg `nsim` replaces `reps`)
