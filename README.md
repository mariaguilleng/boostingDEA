# boostingDEA

This R package, besides implementing the most popular DEA and FDH models, includes two different boosting algorithms for estimating production frontiers: an adaptation of the Gradient Tree Boosting known as EATBoosting and the adaptation of the LS-Boosting algorithm using adapted Multivariate Adaptive Regression Splines (MARS) models as base learners (from now on referred as MARSBoosting). EATBoosting shares similarities with FDH since graphically both generate a step function, while MARSBoosting resembles DEA. However, both algorithms overcome the overfitting problems that characterize standard techniques. Furthermore, in this package, different technical efficiency measures can be calculated. In particular, the input and output-oriented radial measures \citep{banker1984}, the input and output-oriented Russell measures, the Directional Distance Function (DDF), the Weighted Additive Measure (WAM)  and the Slacks-Based Measure (SBM) are included.

## Installation

You can install the released version of the package from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("boostingDEA")
```

And the development version from
[GitHub](https://github.com/itsmeryguillen/boostingDEA) with:

``` r
devtools::install_github("itsmeryguillen/boostingDEA")
```
