<!-- badges: start -->
[![CRAN](http://www.r-pkg.org/badges/version/SSN2)](https://cran.r-project.org/package=SSN2)
[![cran checks](https://badges.cranchecks.info/worst/SSN2.svg)](https://cran.r-project.org/web/checks/check_results_SSN2.html)
[![R-CMD-check](https://github.com/USEPA/SSN2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/USEPA/SSN2/actions/workflows/R-CMD-check.yaml)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![status](https://joss.theoj.org/papers/10.21105/joss.06389/status.svg)](https://joss.theoj.org/papers/10.21105/joss.06389)
<!-- badges: end -->

# SSN2: Spatial Modeling on Stream Networks

`SSN2` is an R package for spatial statistical modeling and prediction on
stream networks, including models based on in-stream distance.
Models are created using moving average constructions. Spatial linear models,
including explanatory variables, can be fit with (restricted) maximum likelihood.
Mapping and other graphical functions are included. It is the successor to the `SSN` R package. See the `SSN2` website for more: [https://usepa.github.io/SSN2/](https://usepa.github.io/SSN2/).

## Citation

If you use `SSN2` in a formal publication or report, please cite it. Citing `SSN2` lets us devote more resources to it in the future. View the `SSN2` citation by running
```r
citation(package = "SSN2")
```

```
#> 
#> To cite SSN2 in publications use:
#> 
#>   Dumelle M, Peterson EE, Ver Hoef JM, Pearse A, Isaak DJ (2024). SSN2:
#>   The next generation of spatial stream network modeling in R. Journal
#>   of Open Source Software, 9(99), 6389,
#>   https://doi.org/10.21105/joss.06389
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {{SSN2}: The next generation of spatial stream network modeling in {R}},
#>     author = {Michael Dumelle and Erin E. Peterson and Jay M. {Ver Hoef} and Alan Pearse and Daniel J. Isaak},
#>     journal = {Journal of Open Source Software},
#>     year = {2024},
#>     volume = {9},
#>     number = {99},
#>     pages = {6389},
#>     doi = {10.21105/joss.06389},
#>     url = {https://doi.org/10.21105/joss.06389},
#>     publisher = {The Open Journal},
#>   }
```

## Statement of Need

Streams provide vital aquatic services that sustain wildlife, provide drinking and irrigation water, and support recreational and cultural activities.  Data are often collected at various locations on a stream network and used to characterize some scientific phenomenon in the stream. Spatial stream network (SSN) models use a spatial statistical modeling framework to describe unique and complex dependencies on a stream network resulting from a branching network structure, directional water flow, and differences in flow volume. SSN models relate a response variable to one or more explanatory variables, a spatially independent error term (i.e., nugget), and up to three spatially dependent error terms: tail-down errors, tail-up errors, and Euclidean errors. Tail-down errors restrict spatial dependence to flow-connected sites (i.e., water flows from an upstream to a downstream site) and incorporate spatial weights (i.e., additive function) to describe the branching network between them. Tail-up errors describe spatial dependence between both flow-connected and flow-unconnected (i.e., sites that share a common downstream junction but not flow) sites, but spatial weights are not required. Euclidean errors describe spatial dependence between sites based on Euclidean distance and are governed by factors not confined to the stream network like regional geology. The `SSN2` **R** package is designed to help users fit SSN models to their stream network data.


## Installation Instructions

Install and load the most recent approved version from CRAN by running
```r
# install the most recent approved version from CRAN
install.packages("SSN2")
# load the most recent approved version from CRAN
library(SSN2)
```

Install and load the most recent version of`SSN2` from GitHub by running
```r
# Installing from GitHub requires you first install the remotes package
install.packages("remotes")

# install the most recent version from GitHub
remotes::install_github("USEPA/SSN2", ref = "main")
# load the most recent development version from GitHub
library(SSN2)
```

Install and load the most recent development version of`SSN2` from GitHub by running
```r
# Installing from GitHub requires you first install the remotes package
install.packages("remotes")

# install the most recent version from GitHub
remotes::install_github("USEPA/SSN2", ref = "develop")
# load the most recent development version from GitHub
library(SSN2)
```

## Contributing to `SSN2`

We encourage users to report bugs and/or contribute to `SSN2`. For more detail on how to do this, please see our contributing guide (`CONTRIBUTING.md`).

## Getting Help

There are several ways to get help with `SSN2`:

1. Open a GitHub issue [link here](https://github.com/USEPA/SSN2/issues).
2. Email the SSN support team (support@spatialstreamnetworks.com or Dumelle.Michael@epa.gov)
3. Post on a support website like Stack Overflow or Cross Validated. 

## Example Usage

Below we provide a brief example showing how to use `SSN2`. For a thorough introduction to the software, see our introductory vignette [linked here](https://usepa.github.io/SSN2/articles/introduction.html). For a list of all functions available in `SSN2`, see our function reference [linked here](https://usepa.github.io/SSN2/reference/index.html).

We load `SSN2`, copy the `.ssn` object that comes with `SSN2` to the temporary directory, and create stream distance matrices used for modeling by running

```r
library(SSN2)
copy_lsn_to_temp()
path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(path, predpts = "pred1km")
ssn_create_distmat(mf04p, predpts = "pred1km", overwrite = TRUE)
```

We fit and summarize an SSN model explaining summer water temperatue (`Summer_mn`) as a function of elevation (`ELEV_DEM`) and precipitation (`AREAWTMAP`) with a exponential, spherical, and Gaussian structures for the tail-up, tail-down, and Euclidean errors, respectively, by running

```r
ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)
summary(ssn_mod)
```

```
#> 
#> Call:
#> ssn_lm(formula = Summer_mn ~ ELEV_DEM + AREAWTMAP, ssn.object = mf04p, 
#>     tailup_type = "exponential", taildown_type = "spherical", 
#>     euclid_type = "gaussian", additive = "afvArea")
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.6393 -2.0646 -0.5952  0.2143  0.7497 
#> 
#> Coefficients (fixed):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 76.195041   7.871574   9.680  < 2e-16 ***
#> ELEV_DEM    -0.026905   0.003646  -7.379  1.6e-13 ***
#> AREAWTMAP   -0.009099   0.004461  -2.040   0.0414 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Pseudo R-squared: 0.6124
#> 
#> Coefficients (covariance):
#>               Effect     Parameter   Estimate
#>   tailup exponential  de (parsill)  3.800e+00
#>   tailup exponential         range  4.194e+06
#>   taildown spherical  de (parsill)  4.480e-01
#>   taildown spherical         range  1.647e+05
#>      euclid gaussian  de (parsill)  1.509e-02
#>      euclid gaussian         range  4.496e+03
#> 
```

We tidy, glance at, and augment (with diagnostics) the fitted model by running

```r
tidy(ssn_mod, conf.int = TRUE)
```

```
#> # A tibble: 3 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept) 76.2       7.87         9.68 0         60.8    91.6     
#> 2 AREAWTMAP   -0.00910   0.00446     -2.04 4.14e- 2  -0.0178 -0.000356
#> 3 ELEV_DEM    -0.0269    0.00365     -7.38 1.60e-13  -0.0341 -0.0198
```

```r
glance(ssn_mod)
```

```
#> # A tibble: 1 × 9
#>       n     p  npar value   AIC  AICc logLik deviance pseudo.r.squared
#>   <int> <dbl> <int> <dbl> <dbl> <dbl>  <dbl>    <dbl>            <dbl>
#> 1    45     3     7  59.3  73.3  76.3  -29.6     41.9            0.612
```

```r
head(augment(ssn_mod))
```

```
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -1515032 ymin: 2529461 xmax: -1512690 ymax: 2531883
#> Projected CRS: USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
#> # A tibble: 6 × 10
#>   Summer_mn ELEV_DEM AREAWTMAP .fitted .resid   .hat  .cooksd .std.resid pid  
#>       <dbl>    <int>     <dbl>   <dbl>  <dbl>  <dbl>    <dbl>      <dbl> <chr>
#> 1     11.4      1977      940.    14.4  -3.07 0.0915 0.0962       -1.78  1    
#> 2     10.7      1984     1087.    12.9  -2.20 0.114  0.00471      -0.352 2    
#> 3     10.4      1993     1087.    12.7  -2.25 0.0372 0.00724      -0.764 3    
#> 4     10.1      2007     1087.    12.3  -2.18 0.0251 0.00153      -0.427 4    
#> 5     10.1      2009     1087.    12.3  -2.13 0.0374 0.000583     -0.216 5    
#> 6      9.81     2012     1109.    12.0  -2.16 0.0602 0.0150       -0.863 6    
#> # ℹ 1 more variable: geometry <POINT [m]>

```

We make predictions at the prediction sites (`pred1km`) by running `predict()` (or `augment()`:

```r
preds <- predict(ssn_mod, newdata = "pred1km", interval = "prediction")
head(preds)
```

```
#>        fit      lwr      upr
#> 1 14.64383 14.27138 15.01627
#> 2 15.00608 14.65017 15.36198
#> 3 14.79235 14.27414 15.31057
#> 4 14.96884 14.45492 15.48276
#> 5 15.15182 14.73770 15.56595
#> 6 15.12783 14.76358 15.49208
```

## Imported Packages

`SSN2` imports the following **R** packages:

* generics: For exporting generic functions.
* graphics: For visualizations (e.g., `plot()`).
* Matrix: For efficient matrix manipulations.
* RSQlite: For various functions that read and write (e.g., `ssn_create_distmat()`).
* sf: For handling spatial data.
* spmodel: For various modeling functions (e.g., `randcov_initial()`) and generic functions (e.g., `loocv()`).
* stats: For various modeling functions (e.g., `confint()`).
* tibble: For creating tibbles as output for various functions (e.g., `tidy()`).
* utils: For various utility functions.
* withr: For path handling while reading and writing.

## Suggested Packages

`SSN2` suggests the following **R** packages:

* ggplot2: For vignette visualizations.
* knitr: For vignette building.
* rmarkdown: For vignette building.
* sp: For making `SSN` objects from the `SSN` **R** package compatible with `SSN2`.
* statmod: For modeling and simulation of inverse Gaussian data.
* testthat: For unit testing.

## License

This project is licensed under the GNU General Public License, [GPL-3](https://cran.r-project.org/web/licenses/GPL-3).

## EPA Disclaimer

The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.

