# MSstatsResponse

**Statistical modeling of drug–protein interactions and IC50 estimation from chemoproteomics data**

MSstatsResponse is an R-based/Bioconductor package for the analysis of mass spectrometry-based chemoproteomics experiments. It is designed to detect drug–protein interactions and estimate IC50 values using isotonic regression models.

The package is compatible with protein-level data from [MSstats](https://github.com/Vitek-Lab/MSstats) and [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT), and supports a range of experimental designs, including dose–response experiments. It is applicable to chemoproteomics workflows based on Data-Independent Acquisition (DIA), Data-Dependent Acquisition (DDA), and targeted Selected Reaction Monitoring (SRM). While primarily developed for chemoproteomics datasets, the methods in MSstatsResponse can be easily extended to other types of datasets that exhibit monotonic trends.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("Vitek-Lab/MSstatsResponse")
```
