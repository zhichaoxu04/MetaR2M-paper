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

  - Meta
  - mvmeta
  - Metafor
  - aod
  - metagen
  - tidyverse
