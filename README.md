# Does sexual history matter? — Reproducibility materials

Companion to the methodological review [**Does sexual history matter?**](https://datepsych.substack.com/p/does-sexual-history-matter) This is an exploratory observational and methodological investigation, not a preregistered or peer-reviewed replication. No participant-level data are distributed here.

## Quick start: regenerate the published figures without microdata

Unzip the archive and work from its root. Use Python 3.12 or 3.13 and install the dependencies (`python -m pip install -r requirements.txt`). Run:

```sh
python forensic/audit.py
python rarticle/update_figures.py
```

These commands regenerate the forensic statistics, ten charts and a standalone HTML article from the archived aggregate results. Figure constants in the figure scripts are rounded transcriptions of the corresponding outputs or published statistics, not new fits. See `FIGURE_SOURCES.md`. The charts are embedded as PNGs in the HTML. The article can be pasted into Substack; individual charts are also available in `article/figures/`. If the editor strips images, insert the PNGs at the corresponding captions. Add the public GitHub URL to the article's Data and reproducibility paragraph after creating the repository; no repository URL has been invented.

## Source data (obtain separately)

- **Date Psychology**: [Who are the men that demand virgins?](https://datepsychology.com/who-are-the-men-that-demand-virgins/) and [“Body Count” And Sexual Double Standards](https://datepsychology.com/body-count-and-sexual-double-standards/). `forensic/audit.py` uses transcribed aggregate numbers and conditional integer-denominator reconstructions. No respondent-level survey data were made available for this audit. Do not interpret reconstructed subgroup counts as verified individual records.
- **NSFG**: obtain the 2015–2017 and 2017–2019 female/male respondent data and Stata dictionary files, and female/male combined 2011–2019 weights data and Stata dictionaries, from [CDC NSFG documentation](https://www.cdc.gov/nchs/nsfg/documentation/combining-data.html). Select the combined `WGT2015_2019` weight. Put the twelve `.dat` and `.dct` files into `nsfg/input/`. The scripts match file names by patterns (downloaded duplicate-name suffixes are accepted): `2015_2017_FemRespData*`, `2015_2017_FemRespSetup*`, `2015_2017_MaleData*`, `2015_2017_MaleSetup*`, the same four patterns for `2017_2019`, `*FemaleWgtData*`, `*FemaleWgtSetup*`, `*MaleWgtData*`, `*MaleWgtSetup*`. Leave only one file per pattern. The public-use files include top-coded measures; never assume unrestricted counts.
- **NLSY97**: obtain the full 1997–2023 public-use archive `nlsy97_all_1997-2023.zip` from [NLSY data access](https://www.nlsinfo.org/content/access-data-investigator/accessing-data-cohorts) (original archive URL: https://nlsinfo.org/cohort-data/nlsy97_all_1997-2023.zip). Place it at the repository root. It contains `nlsy97_all_1997-2023.csv`. `prepare_nlsy.py` extracts the mapped variables, plus sampling design and background variables, into local files. The original analyses used an earlier broader extract; the supplied extractor selects only the variables consumed by the included analysis code. The full newly consolidated extractor has not been end-to-end rerun in this revision; the constituent extraction and ID/frequency checks were previously executed. The full archive is large. Do not upload it, `respondents.csv`, or any generated `_private` files.
- **Thomas et al. (2025)**: [paper](https://doi.org/10.1038/s41598-025-12607-1), [author correction](https://doi.org/10.1038/s41598-025-17764-x), [OSF archive](https://osf.io/b4e8u/). Obtain the archive with `Data/Study 3.csv` and save it as `b4e8u-osfstorage-archive.zip` at the root. Do not redistribute the archive here.
- **Wolfinger & Perry (2023)**: [paper](https://doi.org/10.1016/j.ssresearch.2022.102800), [author manuscript](https://www.dibbleinstitute.org/wp-new/wp-content/uploads/2022/10/Sexual-Resume-and-Marriage-Rates.pdf). The original scripts were requested but were unavailable at the time of revision. The event-history analysis here is an independently specified reconstruction; published and reconstructed analytical intervals differ substantially.

## Run the analyses from data

### NSFG

```sh
cd nsfg
python 01_import.py
python 02_analysis.py
python 03_sensitivity.py
python 04_interactions.py
python 05_menarche.py
python 06_thresholds.py
cd ..
```

`01_import.py` creates a private `nsfg.pkl`. Regression and mean-gap scripts print aggregate output; archived successful runs are in `phase10/nsfg_reexecution.txt`, `phase10/nsfg_interactions.txt` and `nsfg/phase4_menarche_audit_results.txt`. `06_thresholds.py` was executed for this revision and reproduces the reported threshold values to rounding. It uses valid positive counts for each measure separately; virginity uses HADSEX for age 25–49. The 36 strata and four variance units give 108 design degrees of freedom. NSFG marriage associations are retrospective, not prospective or causal.

### Thomas et al.

```sh
python phase10/thomas_robustness.py
```

Within-person matched-profile declines from four to 36 partners, HC3 intervals, country and orientation sensitivities and exploratory equivalence margins. Count conditions were presented in ascending order, so the count effect cannot be separated cleanly from block order. Results are in `phase10/thomas_robustness.csv` and `phase10/thomas_country_estimates.csv`.

### NLSY97

```sh
python prepare_nlsy.py
python phase9/phase9.py
python phase11/audit.py
python phase11/design_audit.py
python phase11/expanded_audit.py
python phase11/joint_ipcw.py
python phase12/reviewer_sensitivity.py
python phase12/reproduce_earlier_core.py
```

If you already have the prior `respondents.csv` extract, use `phase11/extract_design.py` and `phase11/extract_extra.py` instead of `prepare_nlsy.py` to generate the private design and background files from the full archive. Several scripts reuse a prefix of an earlier script to reconstruct precisely the same cohort; maintain these relative paths and do not change the markers. The analysis can be computationally expensive. `phase12/reproduce_earlier_core.py` independently reconstructs the paired attenuation bootstrap and selected event-history specifications; the bootstrap completed and matched the archived values exactly, but the subsequent regressions hit the execution limit in the final validation run. The original aggregate event-history results and the executed Phase 9 timing comparisons are included. A complete fresh run of every script has not been performed in a clean environment.

The primary estimand is the sex-specific standardized first-marriage risk difference by age 40 among eligible unmarried age-20 respondents with a maximum recorded index of 10+ versus 1–4. Initial sampling weights, piecewise exponential hazard models, expanded baseline controls and 117-stratum / 234-PSU linearization are used. The index is not the verified number of distinct lifetime partners; the lower group can have more distinct partners than its maximum. A recent count that is part of the exposure is not included in the primary outcome adjustment because this would require unsupported combinations. Joint censoring-weight sensitivity has no joint confidence interval. The four-category sensitivity includes 7,289 respondents and changes the standardization population. The corroborated-high sensitivity requires numerical reports of at least ten in two distinct interview rounds; 144 high men and 14 high women qualify. Its estimates are imprecise and the pooled fit contains a large coefficient. The restriction changes the estimand and is not a direct reporting-error validation test.

## Archived and executable material

- `forensic/`: executed arithmetic reconstructions and aggregate outputs.
- `nsfg/`: import and analysis code, new threshold estimates, menarche audit.
- `phase9/`: executed reconstruction, paired timing tests, fixed curves, bridge estimates and audit. Conditional hazards are not cumulative attainment. Some bridge adjustments create structural dependence and should remain supplementary.
- `phase10/`: executed Thomas reanalysis and NSFG aggregate logs.
- `phase11/`: executed adversarial audit, survey-design estimates, expanded background adjustment, joint IPCW sensitivity, diagnostics.
- `phase12/`: executed reviewer-requested all-category and repeated-high-report sensitivity, absolute-risk intervals, model audit, independent paired-bootstrap output and the attempted independent event-history rerun code.
- `archived_results/`, `event_history_estimates.csv`, `paired_attenuation.csv`: earlier aggregate outputs preserved for transparency, including less favorable and exploratory results. Exact executable scripts or the full original bootstrap stream for every earlier specification (especially the original Phase 5–8 comprehensive runs, post-2015 models and restricted subgroups) are not all present. They should be regarded as archived outputs, not a claim of fully end-to-end reproduction of every historical result. The headline Phase 11/12 analyses and chart regeneration have executable scripts; the archived event-history results support the Figure 9 intervals.
- `article/`: original figure program, figure amendments and all final charts. Regeneration from aggregate outputs does not require raw data.

## Publication and interpretation


The public-use datasets and original papers remain subject to their own terms. The repository deliberately contains no participant-level data, original copyrighted articles, or external archive copies. Review the files and data-access terms before publishing. `.gitignore` prevents common private intermediates from being committed accidentally.
