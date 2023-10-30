<!-- badges: start -->
[![R-CMD-check](https://github.com/USEPA/SSN2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/USEPA/SSN2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# SSN2: Spatial Modeling on Stream Networks

`SSN2` is an R package for spatial statistical modeling and prediction on
stream networks, including models based on in-stream distance.
Models are created using moving average constructions. Spatial linear models,
including explanatory variables, can be fit with (restricted) maximum likelihood.
Mapping and other graphical functions are included. It is an updated version of the
archived `SSN` R package.

# Installation

Install and load the most recent approved version from CRAN by running
```r
# install the most recent approved version from CRAN
install.packages("SSN2")
# load the most recent approved version from CRAN
library(SSN2)
```

Install and load the most recent development version of`SSN2` from GitHub by running
```r
# Installing from GitHub requires you first install the remotes package
install.packages("remotes")

# install the most recent development version from GitHub
remotes::install_github("USEPA/SSN2", ref = "main")
# load the most recent development version from GitHub
library(SSN2)
```

Install the most recent development version of `SSN2` from GitHub with package vignettes by running
```r
install the most recent development version from GitHub with package vignettes
devtools::install_github("USEPA/SSN2", ref = "main", build_vignettes=TRUE)
```

View the vignette in RStudio by running
```r
vignette("introduction", "SSN2")
```

Further detail regarding SSN2 is contained in the package's documentation manual. 

## Citation

If you use SSN2 in a formal publication or report, please cite it. Citing SSN2 lets us devote more resources to it in the future. View the SSN2 citation by running
```{r}
citation(package = "SSN2")
```

```
#> 
#> To cite SSN2 in publications use:
#> 
#>   Dumelle M, Peterson, E, Ver Hoef JM, Pearse A, Isaak D (2023). SSN2: Spatial Modeling on Stream Networks in R. R package version 0.1.0
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {{SSN2}: Spatial Modeling on Stream Networks in {R}},
#>     author = {Michael Dumelle and Erin Peterson and Jay M. {Ver Hoef} and Alan Pearse and Dan Isaak},
#>     year = {2023},
#>     note = {{R} package version 0.1.0},
#>   }
```

## Package Contributions

We encourage users submit GitHub issues and enhancement requests so we may
continue to improve `SSN2`.

## EPA and Disclaimer

The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.

### License

This project is licensed under the GNU General Public License, [GPL-3](https://cran.r-project.org/web/licenses/GPL-3).
