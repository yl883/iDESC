
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iDESC

<!-- badges: start -->
<!-- badges: end -->

The goal of iDESC is to identify differential expressed gene signatures
between two groups in single cell RNA-seq considering subject effect
with zero-inflated negative binomial mixed model.

## Installation

You can install the released version of iDESC from
[GitHub](https://github.com/yl883/iDESC) with:

``` r
devtools::install_github("yl883/iDESC")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(iDESC)
## basic example code
data(IPF_example)
mat=IPF_example$mat
meta=IPF_example$meta
str(meta)
result=iDESC(mat,meta,subject_var="subject",group_var="disease",zp_group_var="disease")
```
