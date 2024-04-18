# MetaR2M-paper


> [!IMPORTANT]  
> For more details or instructions of this package, please refer to [MetaR2M](https://github.com/zhichaoxu04/MetaR2M). Here we only provide the simulation settings and computational details for the paper.

[MetaR2M](https://github.com/zhichaoxu04/MetaR2M): A novel framework for conducting meta-analyses in high-dimensional settings, specifically for evaluating the R2-based total mediation effect.

We provide the scripts to perform:

- Fixed/Random-effects models simulations with independent/correlated mediators

- Real data analysis with the TOPMed program



## Getting Started

Download and install following required R packages:

- Download [MetaR2M](https://github.com/zhichaoxu04/MetaR2M) package from
  Github using:

<!-- -->

    git clone https://github.com/zhichaoxu04/MetaR2M.git

- Or, install [MetaR2M](https://github.com/zhichaoxu04/MetaR2M) package in
  R directly

  - First, install [devtools](https://devtools.r-lib.org) in R from
    CRAN:

    ``` r
    install.packages("devtools")
    ```

  - Then, install [MetaR2M](https://github.com/zhichaoxu04/MetaR2M) using
    the `install_github` function and load the package:

    ``` r
    devtools::install_github("zhichaoxu04/MetaR2M")
    library(MetaR2M)
    ```

- Make sure that all the required packages have been installed or
  updated. Here are some of the required packages:

  - [Meta](https://cran.r-project.org/web/packages/meta/index.html): User-friendly general package providing standard methods for meta-analysis
  - [mvmeta](https://cran.r-project.org/web/packages/mvmeta/index.html): Collection of functions to perform fixed and random-effects multivariate and univariate meta-analysis and meta-regression.
  - [Metafor](https://www.metafor-project.org/doku.php/metafor): The package consists of a collection of functions that allow the user to calculate various effect size or outcome measures, fit equal-, fixed-, random-, and mixed-effects models to such data, carry out moderator and meta-regression analyses, and create various types of meta-analytical plots.
  - [aod](https://cran.r-project.org/web/packages/aod/index.html): Provides a set of functions to analyse overdispersed counts or proportions.
  - [metagen](http://cran.nexr.com/web/packages/metagen/index.html): Provides methods for making inference in the random effects meta regression model such as point estimates and confidence intervals for the heterogeneity parameter and the regression coefficients vector. 
  - [tidyverse](https://www.tidyverse.org): The tidyverse is an opinionated collection of R packages designed for data science. All packages share an underlying design philosophy, grammar, and data structures.
