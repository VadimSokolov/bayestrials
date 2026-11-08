# ch07_priors.py to R (priors.R)
- events_for_power: dropped, use `ceiling(events_schoenfeld(hr))` (theta_a is `-log(hr)`)
- normal_posterior: dropped, use core `normal_update` (returns list `mean`, `sd`)
- ess_events: `ess_events`
- ess_beta: dropped, core `ess_beta`
- mixture_posterior: `mixture_posterior` (comps is a data.frame with `weight`, `mean`, `sd`; returns list `weights`, `means`, `sds`, `p_pos`)
- cmp_samples(..., rng): `cmp_samples(..., seed)` (returns data.frame `b0`, `slope`)
- prob_success: `prob_success`
