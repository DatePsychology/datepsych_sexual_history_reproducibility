# Figure source map

| Figure | Source of displayed quantities |
|---|---|
| 1 | Hypothetical arithmetic: 104, 27 and 22 no-ideal cases in two distributions of 462; not observed group sizes. |
| 2 | `nsfg/threshold_estimates.csv` (OPPLIFENUM, positive valid counts, age 25–49). Target-sex benchmarking assumes opposite-sex rejection cases. |
| 3 | `phase10/thomas_robustness.csv`, min_pairs=1, specified orientation and covariate/weighting models. |
| 4 | `forensic/welch_reconstruction.csv`; `phase10/nsfg_reexecution.txt` (age 30–49 current-marriage mean differences and direct sex contrast). |
| 5 | `phase10/nsfg_reexecution.txt`, categorical adjusted models 15+ vs 2–4; design-based intervals. |
| 6 | `phase12/reviewer_sensitivity.csv`, analysis=primary_two_groups, estimand=risk. Absolute-risk intervals separately calculated; not borrowed from difference intervals. |
| 6b | Same file, analysis=all_categories, estimand=risk; standardization over all eligible cohort members within sex. |
| 7 | `phase11/expanded_design_estimates.csv`, four index definitions. |
| 8 | `paired_attenuation.csv`; independently reproduced in `phase12/paired_attenuation_independent.csv`. Paired respondent bootstrap, not survey-design intervals. |
| 9 | `event_history_estimates.csv`, common_current and common_lagged, M6_expanded_background, recent-partner coefficients and exponentiated limits; paired differences in `phase9/paired_timing.csv`. |

All figures are original visualizations. Reproduction scripts are `revised_article/make_figures.py` and `update_figures.py`. The original intermediate program generates older versions of certain charts first, then the amendments overwrite them. Running `update_figures.py` is required for the final version.
