# TNBC-Triclosan reproducibility archive

This repository accompanies the manuscript **“A predefined triclosan-related 24-gene score is associated with an EMT-enriched tumor state in triple-negative breast cancer.”**

## Scope

This archive preserves the recovered final/accepted R scripts and provenance files supporting the principal external-validation, direct-exposure, survival, cross-cohort meta-analysis, Hallmark-panorama, and matched-random analyses. It is intended as a transparent analysis archive, not as a claim that every historical intermediate object or the original upstream candidate-collection script has been recovered.

The complete recovered script tree is distributed as `TNBC_Triclosan_repository_ready.zip`. Key provenance files are also exposed directly under `provenance/` for auditability.

## Analysis freeze

The 24-gene panel and downstream inferential analyses are frozen. No gene was replaced, no threshold was altered, no random seed was changed, and no random set was removed based on its result during the recovery/deposition process.

Key final versions include:

- TCGA EMT revalidation: `Script13G_B_v4_COMPACT_FINAL_TCGA_TNBC_EMT_revalidation.R`
- Three-cohort EMT meta-analysis: `Script13G_C_v2_COMPACT_FINAL_cross_cohort_EMT_meta_analysis.R`
- Hallmark preflight: `Script14A_v2_COMPACT_FINAL_Hallmark_panorama_preflight_and_input_lock.R`
- Hallmark cohort analysis/meta-analysis: `Script14B...` + `Script14C...`
- Matched-random analysis: `Script15A...` + `Script15B...` + `Script15C_v2...`

`Script15C v1` is intentionally excluded because it was superseded by v2 after identification of a non-exchangeable GSVA scoring-context/estimand issue.

## Archived candidate resource and 24-gene derivation

The historical raw candidate file contains 251 deduplicated entries under one archived source label (`multi-database+pug`). The file also contains assay-related strings and therefore should **not** be described as 251 standardized human target genes.

Audit against the gene-name field used in the original TCGA analysis showed:

- 251 archived candidate entries
- 66 entries mapped to the analyzed TCGA gene-name space
- 24 mapped entries met the prespecified differential-expression criteria
- those 24 entries exactly match the locked 24-gene panel

See `provenance/` for the raw file and derivation audit.

## Public source datasets

- TCGA-BRCA: NCI Genomic Data Commons
- GSE31519: NCBI Gene Expression Omnibus
- GSE58812: NCBI Gene Expression Omnibus
- GSE95554: NCBI Gene Expression Omnibus
- GSE118389: NCBI Gene Expression Omnibus

Raw public expression data are not redistributed in this repository.

## Script groups inside the packaged archive

- `scripts/01_GSE31519/`: external validation, endpoint curation, survival, and EMT analyses
- `scripts/02_GSE95554_exposure/`: rat exposure design, orthology mapping, and triclosan-vs-oil testing
- `scripts/03_GSE58812/`: preprocessing, score calculation, EMT, and survival analyses
- `scripts/04_TCGA_cross_cohort/`: TCGA revalidation, cross-cohort meta-analysis, Hallmark panorama, and matched-random specificity testing

## Historical paths

The scripts are preserved as executed and may contain historical local Windows paths such as `C:/Users/.../TNBC_Triclosan`. Before rerunning them on another machine, set or replace the project root while keeping the same directory dependencies and input files. The original path strings were not silently rewritten in this archival copy.

## Environment

The final matched-random audit was executed under R 4.6.1 with GSVA 2.6.2. Individual scripts may require additional CRAN/Bioconductor packages documented in their source code.

## Known recovery limitations

1. The original script that assembled the historical 251-entry candidate file was not recovered.
2. The archived candidate file does not preserve row-level CTD/STITCH/PubChem evidence provenance; it contains one historical source label.
3. Some early large R objects, including the original single-cell Seurat object, were not recovered in the historical library.
4. The recovered script archive is therefore described as supporting the principal reported analyses, not as a complete reconstruction of every exploratory historical step.

## Integrity

The packaged archive contains `manifests/repository_manifest_md5.csv`, listing relative paths, byte sizes, and MD5 checksums for the archived files. Script contents are preserved as recovered.

## Repository and archival status

GitHub repository: https://github.com/jjc-jjc39/TNBC-Triclosan

Versioned archive (v1.0.0): https://doi.org/10.5281/zenodo.22638702

DOI: 10.5281/zenodo.22638702
