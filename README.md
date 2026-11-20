# MSstatsResponse

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsResponse.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsResponse/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsResponse/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsResponse/branch/devel)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsResponse is an R/Bioconductor package for semi-parametric modeling of
protein responses in mass spectrometry-based proteomics experiments. Instead of
assuming a fixed parametric curve shape, it fits monotonic curves with isotonic
regression, tests them with F-tests, and estimates confidence intervals by
bootstrapping. It supports two workflows:

- **Dose-response chemoproteomics:** detect drug-protein interactions and
  estimate IC50 values.
- **Protein turnover:** estimate synthesis curves and half-lives from
  stable-isotope labeling time courses (e.g. heavy leucine), with peptide
  quality weights, confidence scores, and turnover categories for each protein.

MSstatsResponse is applicable to DIA, DDA, TMT, and targeted SRM workflows.

MSstatsResponse is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org](http://msstats.org).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsResponse")
```

The development version can be installed directly from this repository:

```r
remotes::install_github("Vitek-Lab/MSstatsResponse")
```

## Quick Start

### Dose-response chemoproteomics

```r
library(MSstatsResponse)

# Example DIA data, already summarized to protein level and normalized
data("DIA_MSstats_Normalized")

# Standardize the data for dose-response modeling
prepared <- MSstatsPrepareDoseResponseFit(
  data = DIA_MSstats_Normalized$ProteinLevelData,
  dose_column = "dose_nM",
  drug_column = "drug",
  protein_column = "Protein",
  log_abundance_column = "LogIntensities",
  transform_nM_to_M = TRUE
)

# Detect drug-protein interactions via isotonic regression + F-test
results <- doseResponseFit(prepared, increasing = FALSE, transform_dose = TRUE)

# Estimate IC50 values with bootstrap confidence intervals
ic50 <- predictIC50(prepared,
                    increasing = FALSE,
                    transform_dose = TRUE,
                    ratio_response = TRUE,
                    bootstrap = TRUE,
                    n_samples = 1000)
```

### Protein turnover

```r
library(MSstatsResponse)

# Simulated heavy-leucine labeling time course (feature-level MSstats format)
data("simulated_turnover_data")

qc <- calculateQCScore(simulated_turnover_data)
ratios <- calculateTurnoverRatios(simulated_turnover_data)
ratios$protein <- ratios$Protein
ratios$drug <- "Synthesis"
ratios$dose <- ratios$TimeVal
ratios$response <- ratios$H_frac

weighted <- calculatePeptideWeights(ratios)
fit <- doseResponseFit(weighted, weights = weighted$weight, increasing = TRUE,
                       transform_dose = FALSE, precalculated_ratios = TRUE)
conf <- calculateConfidence(weighted, fit, qc, simulated_turnover_data)

# Half-life, turnover category, and confidence tier for each protein
final <- classifyTurnoverProteins(weighted, fit, qc, conf)
```

## Input Formats

MSstatsResponse works on data processed with MSstats or MSstatsTMT, not raw
search-tool output:

| Workflow | Upstream function | Field used |
| --- | --- | --- |
| Dose-response | [MSstats](https://github.com/Vitek-Lab/MSstats) `dataProcess()` or [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) `proteinSummarization()` | `$ProteinLevelData` |
| Protein turnover | [MSstats](https://github.com/Vitek-Lab/MSstats) `dataProcess()` on heavy/light labeled data | `$FeatureLevelData` |

For dose-response analysis, `MSstatsPrepareDoseResponseFit()` accepts any data
frame with columns for the dose or concentration, a protein identifier, the
log-abundance response, and, for multi-drug experiments, a drug identifier.
The turnover functions take the feature-level table directly and let you name
the protein, peptide, label, timepoint, and intensity columns. See the
vignettes for details.

## Documentation

- [Dose-response vignette](vignettes/MSstatsResponse.Rmd): full worked example from summarized data to IC50 estimates
- [Protein turnover vignette](vignettes/MSstatsResponse_Protein_Turnover_Workflow.Rmd): turnover ratios, peptide weights, curve fitting, and half-life estimates on a simulated dataset
- [Official website: msstats.org](http://msstats.org)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsResponse)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsResponse/issues).

## References

If you use MSstatsResponse, please cite:

1. Szvetecz, S. et al. (2026). **MSstatsResponse: Semi-parametric statistical
   model enhances detection of drug-protein interactions in chemoproteomics
   experiments.** *Molecular & Cellular Proteomics*, 25(9), 101637.

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsResponse is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
