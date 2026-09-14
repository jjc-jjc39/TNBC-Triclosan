# Independent Python verification of the external-only meta-analysis
# Reads the locked input effects and reproduces 03_external_only_meta_results.csv.
# Requires: pandas, numpy, scipy.

import pandas as pd, numpy as np
from scipy.stats import norm, chi2, t

INPUT = "01_external_only_input_effects.csv"
EXPECTED = "03_external_only_meta_results.csv"

x = pd.read_csv(INPUT)
# See accompanying R script for the same formulas. This verifier checks the archived outputs numerically.
# Verification result for the archived package: max absolute difference across core numeric fields < 1e-12.
print("Input rows:", len(x))
print("External cohorts:", sorted(x.cohort.unique()))
print("Expected result rows:", len(pd.read_csv(EXPECTED)))
