# ch18_cases.py, ctgov_stats.py to R

- `ALPHA, TAU, B_OBF, FRACS` -> dataset `select_design` (list: alpha, tau, obf_bounds, fractions, d_tot, HR, ...)
- `post_prob` -> `select_post_prob`
- `log_bf` -> `select_log_bf`
- `simulate` -> `select_simulate` (args `nsim`, `seed`, `batch`; returns list with `rows` data.frame; calibrates p* by the empirical (1-alpha) quantile of the null path maximum)
- `month_of_events` -> `select_month_of_events` (vectorized in `d`)
- `look_table` -> `select_look_table` (data.frame)
- `ctgov_stats.py` download sections -> `data-raw/ctgov_download.R` (needs jsonlite; untested against the live API); CSV to dataset step -> `data-raw/ctgov.R` producing `ctgov_counts`, `ctgov_bayesian_records`
- `ctgov_stats.py` figures, by-area/phase/sponsor tables, hand-check sample -> dropped: figures are drawn in the book with ggplot2, tables are `table()` on `ctgov_bayesian_records`
