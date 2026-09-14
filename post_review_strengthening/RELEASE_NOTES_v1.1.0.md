# v1.1.0 — post-review strengthening

This release adds two non-reselective, post hoc strengthening components without changing the locked 24-gene panel or any original downstream result:

1. External-only meta-analysis excluding the TCGA derivation/internal-association cohort.
2. Contemporary provenance audit of the 24 locked genes.

The external-only analysis combines GSE31519 and GSE58812 using their accepted cohort effects and empirical Fisher-z variances. The pooled full-EMT effect is rho=0.414688 (normal-theory 95% CI 0.351-0.474). Because only two external cohorts are available, modified Knapp-Hartung inference is explicitly treated as a small-k uncertainty check (P=0.0546), not as a binary validation gate.

The provenance audit is non-reselective. Twenty of 24 locked genes had independently recoverable experimental evidence linked to triclosan at the gene, protein, receptor, or pathway-marker level. TWIST1, ADRB2, ERG, and MITF were not independently recovered and remain in the panel unchanged to preserve prespecification.

The R implementation and an independent Python numerical verifier are both deposited. The document-build environment did not contain Rscript; the Python recomputation reproduced the archived numeric results to machine precision.
