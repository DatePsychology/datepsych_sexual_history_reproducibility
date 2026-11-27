# Phase 11 adversarial audit (NLSY97)

Executed sensitivity analyses of the previously reconstructed age-20 landmark cohorts. Not preregistered, causal, or a direct replication of Wolfinger & Perry. All reported contrasts compare indices >=10 vs 1–4 and predict first marriage through the stated age among eligible unmarried respondents. Index is not a count of verified distinct lifetime partners. Direct-report-only uses nonnegative numerical interval reports and positive numerical lifetime reports through 2005, without skip-derived counts. Regular interval variant uses 9–15-month intervals or lifetime reports; round 1997 count is a 12-month question. The low group may have more distinct partners than its recorded maximum. All new primary estimates use background standardization, not recent-partner adjustment, because recent count is a component of the index. Expanded background adds 1997 religious preference, residence with both biological parents, region, urbanicity, parental religiosity (missing category). These are not full life-course histories.

## Reproduction
Unzip into the directory containing `respondents.csv`, `nlsy97_all_1997-2023.zip`, and the existing `resolved_variables(3).csv` mapping. This archive includes the required Phase 9 source and mapping, but not respondent-level input data. Requires Python, numpy, pandas, scipy, patsy, statsmodels. Run from this directory:

```
python phase11/extract_design.py
python phase11/extract_extra.py
python phase11/audit.py
python phase11/design_audit.py
python phase11/expanded_audit.py
python phase11/joint_ipcw.py
```

The extraction programs stream the very large CSV within the original NLSY97 ZIP; do not publish the generated `_design_private.npz` or `_extra_private.npz` intermediate files. The design scripts validate respondent IDs, strata/PSUs, and the additional baseline variable marginal frequencies against the original codebook. `audit.py` uses the executed Phase 9 reconstruction prefix. `design_audit.py` aggregates model score and empirical covariate-distribution influence contributions over variance PSUs in the full 117-stratum design, including zero contributions outside the analytical domain, uses t intervals with 117 df. This accounts for baseline survey design but not uncertainty in longitudinal weights or initial cohort selection. `expanded_audit.py` repeats design inference with additional baseline controls. `joint_ipcw.py` uses model-based uncensoring probabilities at the start of observed months, capped inverse weights integrated over person-time and applied to events, and standardizes the outcome hazard model over the original weighted analytical cohort. It reports point estimates only (no joint variance estimation); cap 1 reproduces ordinary background standardization. Censoring models use baseline covariates only, not time-varying predictors. All results observational; estimates/intervals are exploratory and nominal. Survey design variables R14897.00/R14898.00; original Round 1 sampling weights.

## Outputs
- `adversarial_estimates.csv`: measurement/overlap sensitivity, respondent-cluster intervals.
- `support_audit.csv`, `raw_reclassification.csv`: aggregate exposure cross-tabs.
- `survey_design_estimates.csv`: design-adjusted original background models.
- `expanded_design_estimates.csv`: expanded baseline controls with design-adjusted intervals.
- `joint_ipcw_estimates.csv`, `joint_ipcw_audit.csv`: joint background/censoring sensitivity and diagnostics.
- `run.log`, `design_run.log`, `expanded_run.log`, `joint_run.log`, `extract.log`, `extra_extract.log`: successful run logs. Initial temporary CSV extraction was rejected by ID check and replaced with binary intermediary.

No respondent-level data in distributed archive. Original inputs remain necessary. Prior analyses should retain original scripts/results and document this post hoc audit as an amendment, not replace their history.
