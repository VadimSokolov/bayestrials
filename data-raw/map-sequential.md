# Mapping: ch01_boundaries.py, ch10_gs.py -> R/sequential.R (and core-sequential.R)

- ch01_boundaries.crossing_prob -> `gs_crossing_prob` (core)
- ch01_boundaries.obf -> `gs_boundary(K, alpha, "obf")` (core)
- ch01_boundaries.pocock -> `gs_boundary(K, alpha, "pocock")` (core)
- ch01_boundaries.naive_inflation -> `gs_naive_inflation(K, z)` (exact, not grid)
- ch10_gs._nodes, _State -> internal (`gs_nodes`, closure in `gs_crossing_two_sided`)
- ch10_gs.crossing -> `gs_crossing_two_sided(t, bounds, lower, drift)` (returns list up, down)
- ch10_gs.type1 -> `gs_power(t, bounds, lower, drift = 0)`; also power at any drift
- ch10_gs.constant_boundary -> `gs_boundary(K, alpha, type, t)` (core, no futility lower bound)
- ch10_gs.spending_boundary -> `gs_spending_boundary` (core)
- ch10_gs.ld_obf, ld_pocock -> dropped: internal to `gs_spending_boundary`
- ch10_gs.power_drift -> `gs_power_drift`
- ch10_gs.expected_information -> `gs_expected_information`
- three tolls table (10-sequential.qmd inline code) -> `gs_tolls`
