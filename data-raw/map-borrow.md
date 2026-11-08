# ch16_borrow.py to R (borrow.R)

- `post_control_components` -> `borrow_post_control` (comps are data frames with w, a, b)
- `prob_superior` (ch16) -> `borrow_prob_superior` (avoids the core `prob_superior`)
- `priors` -> `borrow_priors` (returns a named list of data frames)
- `decision_table` -> `borrow_decision_table`
- `reject_prob` -> `borrow_reject_prob`
- `robust_weight` -> `borrow_robust_weight`
- `ess_robust` -> `borrow_ess_robust`
- `basket_posterior` -> `borrow_basket_posterior` (returns list `means`, `tau_mean`)
- `map_prior_normal` -> `borrow_map_prior` (returns list `samples`, `beta`, `taus`, `pw`; R RNG differs from numpy)
- new: `borrow_power_prior(yh, nh, a0)`
