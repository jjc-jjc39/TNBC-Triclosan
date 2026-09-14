# Independent Python verification of the external-only meta-analysis
# Reads the locked input effects, recomputes all archived result rows, and
# checks agreement with 03_external_only_meta_results.csv.
# Requires: pandas, numpy, scipy.

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2, t

INPUT = "01_external_only_input_effects.csv"
EXPECTED = "03_external_only_meta_results.csv"

x = pd.read_csv(INPUT)
expected = pd.read_csv(EXPECTED)

required = {
    "cohort", "effect_set", "sample_number", "rho", "fisher_z",
    "empirical_fisher_z_variance", "conventional_fisher_z_variance"
}
assert required.issubset(x.columns)
assert len(x) == 4
assert set(x.cohort) == {"GSE31519", "GSE58812"}
assert set(x.effect_set) == {"FULL_EMT", "NO_OVERLAP_EMT"}
assert np.max(np.abs(x.fisher_z.to_numpy() - np.arctanh(x.rho.to_numpy()))) < 1e-12

rows = []

def meta_one(d, variance_col, variance_method):
    z = d.fisher_z.to_numpy(float)
    v = d[variance_col].to_numpy(float)
    k = len(z)
    w = 1.0 / v
    z_fixed = np.sum(w * z) / np.sum(w)
    se_fixed = np.sqrt(1.0 / np.sum(w))
    Q = np.sum(w * (z - z_fixed) ** 2)
    qdf = k - 1
    q_p = chi2.sf(Q, qdf)
    I2 = max(0.0, (Q - qdf) / Q) * 100 if Q > 0 else 0.0
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - qdf) / C) if C > 0 else 0.0
    wr = 1.0 / (v + tau2)
    z_random = np.sum(wr * z) / np.sum(wr)
    se_random = np.sqrt(1.0 / np.sum(wr))
    scale = max(1.0, np.sum(wr * (z - z_random) ** 2) / (k - 1))
    se_mkh = np.sqrt(scale / np.sum(wr))

    def make_row(method, est, se, dist):
        if dist == "NORMAL":
            crit = norm.ppf(0.975)
            p = 2 * norm.sf(abs(est / se))
            df = np.nan
        else:
            df = k - 1
            crit = t.ppf(0.975, df)
            p = 2 * t.sf(abs(est / se), df)
        lo = est - crit * se
        hi = est + crit * se
        return {
            "effect_set": d.effect_set.iloc[0],
            "variance_method": variance_method,
            "pooling_method": method,
            "cohort_number": k,
            "total_sample_number": int(d.sample_number.sum()),
            "pooled_fisher_z": est,
            "standard_error": se,
            "inference_distribution": dist,
            "inference_df": df,
            "pooled_rho": np.tanh(est),
            "confidence_interval_lower_rho": np.tanh(lo),
            "confidence_interval_upper_rho": np.tanh(hi),
            "two_sided_p": p,
            "Cochran_Q": Q,
            "Q_df": qdf,
            "Q_p": q_p,
            "I2_percent": I2,
            "tau2": tau2,
            "modified_KH_residual_scale": scale if method == "RANDOM_DL_MODIFIED_KH" else np.nan,
            "all_external_directions_positive": bool((d.rho > 0).all()),
        }

    return [
        make_row("FIXED_INVERSE_VARIANCE_NORMAL", z_fixed, se_fixed, "NORMAL"),
        make_row("RANDOM_DL_NORMAL", z_random, se_random, "NORMAL"),
        make_row("RANDOM_DL_MODIFIED_KH", z_random, se_mkh, "T"),
    ]

for effect_set, d in x.groupby("effect_set"):
    rows.extend(meta_one(d, "empirical_fisher_z_variance", "EMPIRICAL_BOOTSTRAP_FISHER_Z"))
    rows.extend(meta_one(d, "conventional_fisher_z_variance", "CONVENTIONAL_1_OVER_N_MINUS_3"))

observed = pd.DataFrame(rows)
key_cols = ["effect_set", "variance_method", "pooling_method"]
observed = observed.sort_values(key_cols).reset_index(drop=True)
expected = expected.sort_values(key_cols).reset_index(drop=True)

numeric_cols = [
    "pooled_fisher_z", "standard_error", "pooled_rho",
    "confidence_interval_lower_rho", "confidence_interval_upper_rho",
    "two_sided_p", "Cochran_Q", "Q_p", "I2_percent", "tau2"
]
max_diff = 0.0
for col in numeric_cols:
    diff = np.nanmax(np.abs(observed[col].to_numpy(float) - expected[col].to_numpy(float)))
    max_diff = max(max_diff, float(diff))

assert max_diff < 1e-12, f"Verification failed; max absolute difference={max_diff:.3e}"
print(f"PASS: archived results reproduced; max absolute difference={max_diff:.3e}")
