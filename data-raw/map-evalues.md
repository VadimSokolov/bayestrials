# Mapping: ch11_sim.py -> R/evalues.R

- mix_eprocess -> `ev_mixture`
- log_mix_eprocess -> `ev_log_mixture`
- walks (generator) -> `ev_walks` (returns one matrix; callers chunk)
- first_cross -> `ev_first_cross`
- naive_false_alarm -> `ev_naive_false_alarm` (returns list naive, evalue)
- select_like -> `ev_select_like` (list cb, co, nmax, null, alt)
- betting_paths -> `ev_betting_paths`
