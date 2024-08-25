This repository contains all simulation and example code for the paper "Generalized Confidence Bounds for Multiple Testing, with Application to Adaptive False Discovery Rate". To use the code, you must install this package by

```r
devtools::install_github("Jiefei-Wang/Confidence_envelope_paper")
```

Here is a list of the required packages for the code:
```r
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("cp4p")

## for bioconductor packages
install.packages("BiocManager")
BiocManager::install("BiocParallel")

## for github packages
install.packages("devtools")
devtools::install_github('jtleek/tidypvals')
```

Below is the files in this repository:

- `simulation1.r`: Estimate the coverage probability
- `simulation2.r`: Estimate Expectation and MSE
- `simulation3.r`: Estimate the FDR
- `example.r`: Code for the paper example

