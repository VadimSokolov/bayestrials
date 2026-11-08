# Mapping: code/ch13_predictive.py to R/predictive.R

- `ZA` -> dropped: use `qnorm(0.975)` (default `za` of `cond_power`, `pp_normal`)
- `cond_power` -> `cond_power`
- `pred_prob_normal` -> `pp_normal`
- `pp_one_arm` -> `pp_one_arm`
- `lee_liu_table` -> `lee_liu_table` (returns data.frame `n, s_fut, s_eff` instead of a dict)
- `lee_liu_oc` -> `lee_liu_oc` (returns named vector `P_succ, EN`)
- `_grid`, `pr_superior` -> internal `.gl_grid`, `.pr_superior_grid` (use `prob_superior()` for general use)
- `pp_two_arm` -> `pp_two_arm`
- `backward_induction` -> `backward_induction` (list `stop`, `action`; element k+1 is look k)
- `forward_oc` -> `forward_oc` (takes the `backward_induction` result, returns `c(P, EN)`)
- `fixed_oc` -> `fixed_oc`
