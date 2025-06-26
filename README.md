# MSstatsResponse

**Statistical Modeling of Drug–Protein Interactions and IC50 Estimation from Dose–Response Proteomics Data**

MSstatsResponse is an R package for analyzing mass spectrometry-based dose-response proteomics experiments. It enables robust detection of drug–protein interactions and estimation of IC50 values using isotonic regression.

The package is compatible with protein-level data from [MSstats](https://github.com/Vitek-Lab/MSstats) and [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT), and is designed to support a wide range of experimental designs.

## Features

- Drug–protein interaction detection using isotonic regression and F-tests  
- IC50 estimation via monotonic model fitting and linear interpolation  
- Bootstrap confidence intervals for IC50 estimates  
- Publication-ready plots of dose–response curves  
- Supports both log2 and ratio scale analysis, and flexible monotonicity constraints

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("xx/MSstatsResponse")
```
