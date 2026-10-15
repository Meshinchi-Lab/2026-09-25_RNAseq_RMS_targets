## Rhabdomyosarcoma Therapeutic Targets — `rms_gene_explorer`

A Shiny app for data mining gene expression from rhabdomyosarcoma tumors and cell lines.
Select one gene; see its distribution across every cohort assembled for the RMS CAR-T
target-discovery project.

### Datasets

| Dataset | n samples | Assay |
|---|---|---|
| PCGP (St. Jude) | 169 | pTPM |
| dGTEx (pediatric normal) | 332 | pTPM |
| Adult GTEx v11 (all 30 tissues) | 19,616 | pTPM |
| GSE108022 (Shern 2014) | 106 | pTPM |
| DepMap Public 26Q1 (14 RMS-FP, 11 RMS-FN, 1,694 other) | 1,719 | log2(TPM+1) |
| Wachtel 2023 RMS surfaceome | 15 | iTop3 abundance |

Source objects are read in place from `2025-12-04_RNAseq_RMS_infant/results/`; set
`source_dir` in `config.yml`. Nothing is copied into this repo.

### Build the store

```bash
Rscript bin/01_build_store.R        # melt sources to staging/ (~2 min)
Rscript bin/02_repartition_store.R  # repartition by gene bucket
Rscript bin/03_build_indexes.R      # samples.parquet + gene_index.parquet
Rscript bin/04_validate_store.R     # reconcile against sources
```

Stage 1 accepts dataset names to build incrementally, e.g.
`Rscript bin/01_build_store.R GSE108022 Wachtel_protein`.

### Run

Run locally with bare metal installation of R and package dependencies for the app installed. 

```bash
Rscript -e 'shiny::runApp("shiny/rms_gene_explorer", port = 8080)'
```

Or use the docker image once it is built or pulled from the container registry. 

```bash
docker run --rm \
    -v "$PWD":"$PWD" \
    -v "$PWD/store":"$PWD/store":ro \
    -w "$PWD" \
    -p 127.0.0.1:3838:3838 \
    ghcr.io/meshinchi-lab/rms-gene-explorer:4.5.3
```

### Tests

```bash
Rscript tests/testthat.R
```

### Panels

Three tabs share one gene selection:

1. **RMS vs healthy muscle** — PCGP, GSE108022 and DepMap against skeletal muscle only.
   dGTEx and adult GTEx contribute their 48 and 818 Muscle samples; every other tissue is
   excluded here so the normal median is not dragged toward whichever tissue dominates
   the cohort.
2. **Normal Tissues** — the on-target/off-tumour view, with a scope toggle. *Critical
   tissues + muscle* (default) shows the CAR-T safety-critical set: 9 tissues in adult
   GTEx, 8 in dGTEx, which has no whole-blood tissue. *All tissues* shows all 24 / 30.
   This is why all 30 GTEx tissues are stored rather than a critical subset.
3. **Protein (Wachtel)** — iTop3 abundance on its own axis, coloured by model class
   (FP cell line / FN cell line / PDX / control) with fusion status on hover.

DepMap ships values that are already `log2(TPM+1)`. The scale toggle is assay-aware, so
DepMap is never logged twice and is back-transformed on the linear setting. 

### Data Sources

The expression store is narrow long-format Parquet (`gene_symbol`, `dataset`,
`sample_id`, `value`) hive-partitioned on the first letter of the gene symbol, so a
single-gene query opens one directory instead of scanning ~411M rows. Sample metadata
lives in a separate ~22k-row table held in memory and joined at query time — not
repeating group and tissue strings across every row is what keeps the store near 1.5 GB.

`R/plot_gene_distribution.R` is called by both this app and
`2026-08-25_RMS_multi_cohort_expression_figures.qmd` found at the [companion repository](https://github.com/Meshinchi-Lab/2025-12-04_RNAseq_RMS_infant). 

Adding a dataset (R2 microarray, HPA mass-spec, CCLE RPPA) is supported by `register_loader()` in `R/loaders.R`. 

### Data governance

PCGP is controlled-access under dbGaP **phs000720** and the OPBG cohort is pre-publication. Run only on the analysis host or behind institutional SSO. **Do not deploy to shinyapps.io and do not package with shinylive** — shinylive ships the entire store into the viewer's browser. `store/`, `staging/`, `meta_staging/`, and `data/` are gitignored; no sample-level values belong in this repo.


Author: Jenny L. Smith
ORCID: [0000-0003-0402-2779](https://orcid.org/0000-0003-0402-2779)
