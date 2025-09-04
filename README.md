# MSstatsResponse

**Statistical modeling of drug–protein interactions and IC50 estimation from chemoproteomics data**

MSstatsResponse is an R package for analyzing mass spectrometry-based chemoproteomics proteomics experiments. It uses isotonic regression to detect drug–protein interactions and estimate IC50 values.

The package is compatible with protein-level data from [MSstats](https://github.com/Vitek-Lab/MSstats) and [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT), and is designed to support a wide range of experimental designs.

## Features

- Drug–protein interaction detection using isotonic regression and F-tests  
- IC50 estimation via monotonic model fitting and linear interpolation  
- Bootstrap confidence intervals for IC50 estimates  
- Dose-response curve plots 
- Supports both log2 and ratio scale analysis, and flexible monotonicity constraints
- Can be extended to use for analysis of any datasets that exhibit monotonicity trends

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("Vitek-Lab/MSstatsResponse")
```
