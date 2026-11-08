# bayestrials

Companion R package to the book *Bayesian Adaptive Clinical Trials* by Nicholas G. Polson and Vadim Sokolov. It implements every computation in the book, and the book's code chunks are short calls to it.

## Installation

```r
# install.packages("remotes")
remotes::install_github("VadimSokolov/bayestrials")
```

## What it does

| Family | Functions (examples) | Book chapter |
|---|---|---|
| Posterior computation | `beta_post`, `prob_superior`, `normal_update`, `ess_beta`, `ess_normal`, `hpd_beta`, `dbetabinom` | Prior, likelihood, posterior |
| Testing and randomization basics | `z_test_binary`, `wald_ci`, `perm_test`, `assign_block`, `assign_minimization`, `arnold_stochasticity` | Testing; Randomization |
| Evidence | `bf10_normal`, `min_bf`, `false_positive_risk`, `calibrate_c_flat` | Weighing evidence |
| Power and sample size | `n_fixed`, `events_schoenfeld`, `assurance`, `n_md`, `events_md` | Error rates, power, sample size |
| Multiplicity and shrinkage | `mult_holm`, `mult_graph`, `shrink_eb_fit`, `shrink_post_prob` | Multiplicity |
| Group sequential designs | `gs_boundary`, `gs_spending_boundary`, `gs_crossing_prob`, `gs_power`, `gs_tolls` | Interim looks |
| E-values | `ev_mixture`, `ev_first_cross`, `ev_betting_paths` | Anytime-valid inference |
| Operating characteristics | `oc_simulate`, `oc_calibrate`, `oc_summarize` | Operating characteristics by simulation |
| Predictive probability | `pp_two_arm`, `cond_power`, `lee_liu_table`, `backward_induction` | Predictive probability |
| Adaptive randomization | `rar_simulate`, `rar_gittins_index`, `rar_ecmo_urn` | Response-adaptive randomization |
| Dose finding | `dose_three_plus_three`, `dose_crm`, `dose_boin` | Dose finding |
| Borrowing | `borrow_power_prior`, `borrow_map_prior`, `borrow_ess_robust` | Borrowing information |
| Platform trials | `platform_shared_control_table`, `platform_time_machine` | Master protocols and platforms |
| Case studies | `select_simulate`, `select_look_table`; datasets `ulcer`, `select_design`, `ctgov_counts` | Case studies |

A sponsor's calibration workflow in four calls:

```r
library(bayestrials)
looks <- c(40, 80, 120, 160)
null  <- oc_simulate(0.30, 0.30, looks, nsim = 1e5, seed = 11)  # simulate under the null
c0    <- oc_calibrate(null, alpha = 0.025)                    # posterior threshold for 2.5% type I error
alt   <- oc_simulate(0.30, 0.45, looks, nsim = 4e4, seed = 12)
oc_summarize(alt, c0)                                         # power, expected sample size, early stopping
```

Methods are checked against published values (for example the O'Brien-Fleming boundaries 4.049, 2.863, 2.337, 2.024 and the Lan-DeMets spending boundaries of gsDesign). For production design work see also rpact, gsDesign, RBesT and BOIN on CRAN.

## Interactive companion

https://nytel.ai/testing/ (including a browser-based Bayesian trial simulator).

## License

MIT
