# Mapping: code/ch12_sim.py to R/simulate.R

- `post_table` -> `oc_post_table`
- `simulate` -> `oc_simulate` (adds `method = "normal"` option; `nrep` -> `nsim`)
- `with_futility` -> `oc_with_futility`
- `predictive_table` -> `oc_predictive_table`
- `predictive_matrix` -> `oc_predictive_matrix`
- `max_pp` -> `oc_max_pp`
- `calibrate` -> `oc_calibrate` (returns the calibrated c)
- `summarize` -> `oc_summarize` (unused `fut` argument dropped)
