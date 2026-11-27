NSFG 2015–2019 exploratory marriage-outcome analysis

Requirements: Python 3, pandas, numpy, scipy, patsy.
Extract the supplied CDC.zip directly into input/ beside these scripts (original filenames). From this directory run, in order:
python 01_import.py
python 02_analysis.py
python 03_sensitivity.py
python 04_interactions.py

Scripts print aggregate results; nsfg.pkl is a local intermediate containing respondent data and is not included. No microdata distributed. The survey variance estimator uses stratum-centered PSU score totals including zero contributions outside each domain, with 108 design degrees of freedom. Primary outcome is LIFPRTNR (50+ coded 50). Analyses are exploratory cross-sectional associations, not causal estimates. PTSB4MAR analysis is limited to ever-married women and includes the first spouse if premarital intercourse occurred.
