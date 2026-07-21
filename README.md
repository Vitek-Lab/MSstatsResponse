# MSstatsResponse

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsResponse.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsResponse/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsResponse/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsResponse/branch/devel)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsResponse is an R/Bioconductor package for semi-parametric protein
response modeling of mass spectrometry-based proteomics experiments. It is
designed for dose-response chemoproteomics — detecting drug-protein interactions
and estimating IC50 values — using isotonic regression, F-tests, and
bootstrap-based confidence intervals rather than assuming a fixed parametric
curve shape. The same framework applies to protein turnover kinetics
(synthesis/degradation curves and half-life estimation) and, more generally, to
any assay whose response follows a monotonic trend with a dose, concentration,
time, or temperature. MSstatsResponse is applicable to DIA, DDA, TMT, and targeted
SRM workflows.

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

## Input Formats

MSstatsResponse works on **summarized protein-level data**, not raw search-tool
output. It is compatible with the output of the standard MSstats and MSstatsTMT
summarization steps:

| Upstream package | Summarization function | Field used |
| --- | --- | --- |
| [MSstats](https://github.com/Vitek-Lab/MSstats) | `dataProcess()` | `$ProteinLevelData` |
| [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) | `proteinSummarization()` | `$ProteinLevelData` |

`MSstatsPrepareDoseResponseFit()` accepts any data frame with columns for the
dose/concentration/time (x-axis), a protein/gene/peptide identifier, the
log-abundance response (y-axis), and — for multi-drug experiments — a drug
identifier. This makes the workflow extensible to other assays with a monotonic
response, including protein turnover kinetics. See the vignette for details on
each required column.

## Documentation

- [MSstatsResponse vignette](vignettes/MSstatsResponse.Rmd) — full worked example from summarized data to IC50 estimates
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

1. Szvetecz S, Kohler D, Vitek O. **MSstatsResponse: Semi-parametric statistical
   model enhances detection of drug-protein interactions in chemoproteomics
   experiments.** *bioRxiv*. 2026.
   [DOI: 10.64898/2026.03.09.710598](https://doi.org/10.64898/2026.03.09.710598)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsResponse is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
