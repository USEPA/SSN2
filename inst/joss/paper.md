---
title: 'SSN2: The next generation of spatial stream network modeling in R'
tags:
  - R
  - Spatial stream network
  - Spatial correlation
  - Spatial prediction
  - Geostatistics
  - Generalized linear model
  - Torgegram
authors:
  - name: Michael Dumelle
    orcid: 0000-0002-3393-5529
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Erin E. Peterson
    affiliation: 2
  - name: Jay M. Ver Hoef
    affiliation: 3
  - name: Alan Pearse
    affiliation: 4
  - name: Daniel J. Isaak
    affiliation: 5
affiliations:
 - name: Office of Research and Development, United States Environmental Protection Agency (USEPA)
   index: 1
 - name: EP Consulting and Centre for Data Science, Queensland University of Technology, Brisbane Australia 4000 
   index: 2
 - name: NMFS Alaska Fisheries Science Center, United States National Oceanic and Atmospheric Administration (NOAA)
   index: 3
 - name: NIASRA, School of Mathematics and Applied Statistics, University of Wollongong
   index: 4
 - name: Rocky Mountain Research Station, United States Forest Service (USFS)
   index: 5
citation_author: Dumelle et. al.
date: 07 July 2024
year: 2024
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
editor_options: 
  chunk_output_type: console
---

# Summary

The `SSN2` **R** package provides tools for spatial statistical modeling, parameter estimation, and prediction on stream (river) networks. `SSN2` is the successor to the `SSN` **R** package [@ver2014ssn], which was archived alongside broader changes in the **R**-spatial ecosystem [@nowosad2023] that included 1) the retirement of `rgdal` [@bivand2021rgdal], `rgeos` [@bivand2020rgeos], and `maptools` [@bivand2021maptools] and 2) the lack of active development of `sp` [@bivand2013applied]. `SSN2` maintains compatibility with the input data file structures used by the `SSN` **R** package but leverages modern **R**-spatial tools like `sf` [@pebesma2018sf]. `SSN2` also provides many useful features that were not available in the `SSN` **R** package, including new modeling and helper functions, enhanced fitting algorithms, and simplified syntax consistent with other **R** generic functions.

# Statement of Need

Streams provide vital aquatic services that sustain wildlife, provide drinking and irrigation water, and support recreational and cultural activities.  Data are often collected at various locations on a stream network and used to characterize spatial patterns in stream phenomena. For example, a manager may need to know how the amount of a hazardous chemical changes throughout a stream network to inform mitigation efforts. Comprehensive formulations of spatial stream network (SSN) models are provided by @ver2010moving, @peterson2010mixed, and @ver2014ssn. The `SSN2` **R** package is designed to help users fit SSN models to their stream network data.

SSN models use a spatial statistical modeling framework [e.g., @cressie1993statistics] to describe unique and complex dependencies on a stream network resulting from a branching network structure, directional water flow, and differences in flow volume. These SSN models relate a continuous or discrete response variable to one or more explanatory variables, a spatially independent random error term (i.e., nugget), and up to three spatially dependent random error terms: tail-up random errors, tail-down random errors, and Euclidean random errors. Tail-up random errors restrict spatial dependence to flow-connected sites (i.e., water flows from an upstream to a downstream site) and incorporate spatial weights through an additive function to describe the branching network between sites. Tail-down random errors describe spatial dependence between both flow-connected and flow-unconnected sites (i.e., sites that share a common downstream junction but not flow), but spatial weights are not required. Euclidean random errors describe spatial dependence between sites based on straight-line distance and are governed by factors not confined to the stream network, such as regional geology. The variances and the length-scales of spatial dependence in the tail-up, tail-down, and Euclidean random errors are controlled by separate variance (i.e., partial sill) and range parameters, respectively. In this paper, we show how to use the `SSN2` **R** package to fit SSN models, inspect SSN models, and use SSN models to make predictions at unobserved locations on a stream network. 

# Package Overview

The streams, observation, and prediction datasets must be pre-processed prior to fitting SSN models and making predictions at unobserved locations using `SSN2`. Previously, the STARS toolset for ArcGIS Desktop versions 9.3x - 10.8x [@peterson2014stars] or the `openSTARS` **R** package [@kattwinkel2020preparing] were used to generate spatial information required for model fitting and prediction. However, both software packages have recently been retired and are replaced by the `SSNbler` **R** package [@peterson2024SSNbler], which is a new, **R**-based version of the STARS tools. `SSNbler` is currently available on GitHub, will soon be available on CRAN, and contains several useful resources that guide users through these pre-processing steps. Pre-processing using either `SSNbler`, STARS, or `openSTARS` ends with the creatino of a `.ssn` folder, which is non-proprietary. Files residing in the `.ssn` folder are read into R using `ssn_import()` from `SSN2` and palced into a list structure called an SSN object, which contains all the spatial, topological, and attribute information needed to leverage the modeling tools in `SSN2`.

`SSN2` is first installed from CRAN:

``` r
install.packages("SSN2")
```

Then, `SSN2` is loaded into an **R** session:

``` r
library(SSN2)
```

The `SSN2` package comes with an example `.ssn` folder called `MiddleFork04.ssn` that represents water temperatures recorded from a stream network in the Middle Fork of the Salmon River in Idaho, USA during 2004. 

Several functions in `SSN2` for reading and writing data (which we use shortly) directly manipulate the `.ssn` folder. To avoid directly manipulating the `MiddleFork04.ssn` data installed alongside `SSN2`, `MiddleFork04.ssn` is instead copied into a temporary directory and the relevant path to this directory stored:

``` r
copy_lsn_to_temp()
path <- file.path(tempdir(), "MiddleFork04.ssn")
```

The `copy_lsn_to_temp()` function is only used when working with `MiddleFork04.ssn` and generally, `path` should indicate a permanent directory on your computer that points towards your `.ssn` object. After specifying `path`, the stream reaches, observed sites, and prediction sites (`pred1km`) are imported and then visualized (Figure$~$\ref{fig:steam-network}):

``` r
mf04p <- ssn_import(path, predpts = "pred1km")
```


``` r
library(ggplot2)
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$preds$pred1km, pch = 17, color = "blue") +
  geom_sf(data = mf04p$obs, color = "brown", size = 2) +
  theme_bw()
```

\begin{figure}

{\centering \includegraphics[width=0.85\linewidth]{paper_files/figure-latex/steam-network-1} 

}

\caption{Middle Fork 2004 stream networks. Observed sites are represented by brown, closed circles. Prediction sites are represented by blue, closed triangles.}\label{fig:steam-network}
\end{figure}

Prior to statistical modeling, hydrologic distance matrices are created: [@ver2010moving]:

``` r
ssn_create_distmat(mf04p, predpts = "pred1km", overwrite = TRUE)
```

Of particular interest here is summer mean stream temperature (`Summer_mn`) in degrees Celsius, 
which we will model as a function of elevation (`ELEV_DEM`) and watershed-averaged precipitation (`AREAWTMAP`) with exponential, spherical, and Gaussian structures for the tail-up, tail-down, and Euclidean errors, respectively. Using `ssn_lm()`, the model is fit:

``` r
ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)
```

The `additive` argument represents an "additive function value (AFV)" variable that captures branching in the stream network and is required when modeling the tailup covariance. Cumulative watershed area is commonly used to derive the additive function value (here, `afvArea` represents cumulative watershed area), but other variables like flow can be used (if every line feature in the `edges` dataset contains a non-null value). @ver2010moving provide further details regarding additive function values.

The `ssn_lm()` function is designed to be similar in syntax and structure to the `lm()` function in base **R** for fitting nonspatial linear models. Additionally, `SSN2` accommodates various S3 methods for commonly-used **R** generic functions that operate on model objects. For example, the generic function `summary()` is used to summarize the fitted model:

``` r
summary(ssn_mod)
```

```
## 
## Call:
## ssn_lm(formula = Summer_mn ~ ELEV_DEM + AREAWTMAP, ssn.object = mf04p, 
##     tailup_type = "exponential", taildown_type = "spherical", 
##     euclid_type = "gaussian", additive = "afvArea")
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -3.6393 -2.0646 -0.5952  0.2143  0.7497 
## 
## Coefficients (fixed):
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept) 76.195041   7.871574   9.680  < 2e-16 ***
## ELEV_DEM    -0.026905   0.003646  -7.379  1.6e-13 ***
## AREAWTMAP   -0.009099   0.004461  -2.040   0.0414 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Pseudo R-squared: 0.6124
## 
## Coefficients (covariance):
##               Effect     Parameter   Estimate
##   tailup exponential  de (parsill)  3.800e+00
##   tailup exponential         range  4.194e+06
##   taildown spherical  de (parsill)  4.480e-01
##   taildown spherical         range  1.647e+05
##      euclid gaussian  de (parsill)  1.509e-02
##      euclid gaussian         range  4.496e+03
##               nugget        nugget  2.087e-02
```

`SSN2` methods for the `tidy()`, `glance()`, and `augment()` generic functions from the `broom` **R** package [@robinson2021broom] are used to inspect the fitted model and provide diagnostics:

``` r
tidy(ssn_mod, conf.int = TRUE)
```

```
## # A tibble: 3 x 7
##   term        estimate std.error statistic  p.value conf.low conf.high
##   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
## 1 (Intercept) 76.2       7.87         9.68 0         60.8    91.6     
## 2 AREAWTMAP   -0.00910   0.00446     -2.04 4.14e- 2  -0.0178 -0.000356
## 3 ELEV_DEM    -0.0269    0.00365     -7.38 1.60e-13  -0.0341 -0.0198
```

``` r
glance(ssn_mod)
```

```
## # A tibble: 1 x 9
##       n     p  npar value   AIC  AICc logLik deviance pseudo.r.squared
##   <int> <dbl> <int> <dbl> <dbl> <dbl>  <dbl>    <dbl>            <dbl>
## 1    45     3     7  59.3  73.3  76.3  -29.6     41.9            0.612
```

``` r
aug_mod <- augment(ssn_mod)
subset(aug_mod, select = c(Summer_mn, .fitted, .resid, .hat, .cooksd))
```

```
## Simple feature collection with 45 features and 5 fields
## Geometry type: POINT
## Dimension:     XY
## Bounding box:  xmin: -1530805 ymin: 2527111 xmax: -1503079 ymax: 2537823
## Projected CRS: USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
## # A tibble: 45 x 6
##    Summer_mn .fitted .resid   .hat  .cooksd           geometry
##        <dbl>   <dbl>  <dbl>  <dbl>    <dbl>        <POINT [m]>
##  1     11.4     14.4  -3.07 0.0915 0.0962   (-1512690 2531883)
##  2     10.7     12.9  -2.20 0.114  0.00471  (-1512852 2531295)
##  3     10.4     12.7  -2.25 0.0372 0.00724  (-1513400 2530706)
##  4     10.1     12.3  -2.18 0.0251 0.00153  (-1514027 2530147)
##  5     10.1     12.3  -2.13 0.0374 0.000583 (-1514309 2529902)
##  6      9.81    12.0  -2.16 0.0602 0.0150   (-1515032 2529461)
##  7      9.76    11.6  -1.85 0.0736 0.00739  (-1515513 2528810)
##  8      9.77    11.6  -1.84 0.0648 0.00687  (-1515588 2528592)
##  9      9.53    11.4  -1.87 0.112  0.00152  (-1516440 2527899)
## 10     12.6     14.9  -2.28 0.0498 0.00964  (-1512464 2531195)
## # i 35 more rows
```

Specific generic helper functions (e.g., `coef()`, `AIC()`, `residuals()`) can be used to obtain the same quantities returned by `tidy()`, `glance()`, and `augment()`:

``` r
coef(ssn_mod)
```

```
## (Intercept)    ELEV_DEM   AREAWTMAP 
## 76.19504087 -0.02690478 -0.00909941
```

``` r
AIC(ssn_mod)
```

```
## [1] 73.2623
```

``` r
head(residuals(ssn_mod))
```

```
##         1         2         3         4         5         6 
## -3.066413 -2.204147 -2.252004 -2.175337 -2.131527 -2.162417
```

Spatial prediction (i.e., Kriging) at the unobserved sites is performed using the generic functions `predict()` or `augment()`:

``` r
aug_pred <- augment(ssn_mod, newdata = "pred1km", interval = "prediction")
subset(aug_pred, select = c(.fitted, .lower, .upper))
```

```
## Simple feature collection with 175 features and 3 fields
## Geometry type: POINT
## Dimension:     XY
## Bounding box:  xmin: -1530631 ymin: 2521707 xmax: -1500020 ymax: 2540253
## Projected CRS: USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
## # A tibble: 175 x 4
##    .fitted .lower .upper           geometry
##      <dbl>  <dbl>  <dbl>        <POINT [m]>
##  1    14.6   14.3   15.0 (-1520657 2536657)
##  2    15.0   14.7   15.4 (-1519866 2536812)
##  3    14.8   14.3   15.3 (-1521823 2536911)
##  4    15.0   14.5   15.5 (-1523183 2537256)
##  5    15.2   14.7   15.6 (-1523860 2537452)
##  6    15.1   14.8   15.5 (-1525443 2537698)
##  7    15.1   14.7   15.5 (-1526397 2537254)
##  8    15.0   14.6   15.4 (-1527436 2536803)
##  9    14.9   14.6   15.3 (-1529043 2536449)
## 10    14.9   14.5   15.3 (-1529689 2537313)
## # i 165 more rows
```

Here, `.fitted` are the predictions, `.lower` are the lower bounds of 95% prediction intervals, and `.upper` are the upper bounds of 95% prediction intervals. Utilizing `augment()` makes the prediction output straightforward to visualize:


``` r
ggplot() +
    geom_sf(data = mf04p$edges) +
    geom_sf(data = aug_pred, aes(color = .fitted), size = 2) +
    scale_color_viridis_c(name = "Pred.", option = "H") +
    theme_bw()
```

\begin{figure}

{\centering \includegraphics[width=0.85\linewidth]{paper_files/figure-latex/steam-preds-1} 

}

\caption{Predicted Middle Fork 2004 mean summer temperatures (Celsius) spaced one kilometer apart. }\label{fig:steam-preds}
\end{figure}

Generalized spatial linear models for binary, count, proportion, and skewed data are available via the `ssn_glm()` function. `ssn_lm()` and `ssn_glm()` also accommodate several advanced features, which include nonspatial random effects as in `lme4` and `nlme`; see @bates2015lme4 and @pinheiro2006mixed, respectively and Euclidean anisotropy [@zimmerman2024spatial], among others. In addition to modeling, simulating data on a stream network is performed via `ssn_simulate()`. 

# Discussion

SSN models are valuable tools for statistical analysis of data collected on stream networks and help improve inference about vital stream ecosystems. These models have been employed (using `SSN`) to better understand and manage water quality [@scown2017improving; @mcmanus2020variation], ecosystem metabolism [@rodriguez2019estimating], and climate change impacts on freshwater ecosystems [@ruesch2012projected; @isaak2017norwest], as well as generate aquatic population estimates [@isaak2017scalable], inform conservation planning [@rodriguez2019spatial; @sharma2021dendritic], and assess restoration activities [@fuller2022riparian], among other applications. The breadth and applicability of SSN models are further enhanced by data aggregation tools like the National Hydrography Dataset [@mckay2012nhdplus], National Stream Internet Project [@nagel2015national] and StreamCat [@hill2016stream].

There are several spatial modeling packages in **R**, including `geoR` [@ribiero2022geoR], `gstat` [@pebesma2004gstat], `FRK` [@sainsbury2024modeling], `fields` [@nychka2021fields], `R-INLA` [@lindgren2015bayesian], and `spmodel` [@dumelle2023spmodel], among others. However, these aforementioned spatial modeling packages do not account for the unique spatial relationships found in data collected on stream networks. The `rtop` [@skoien2014rtop], `VAST` [@charsley2023catchment], and `SSN2` **R** packages can be used to describe spatial stream network data in **R**, but `SSN2` is unique. It not only provides representations of stream network data in **R** but also provides an extensive suite of functions for model fitting, diagnostics, and spatial integration that integrate with the popular "tidy" framework [@wickham2019welcome; @kuhn2022tidy]. To learn more about `SSN2`, visit the CRAN webpage at [https://CRAN.R-project.org/package=SSN2](https://CRAN.R-project.org/package=SSN2).

# Acknowledgements

Figures were created using `ggplot2` [@wickham2016ggplot2] and the `viridis` color palettes [@garnier2024viridis].

<!-- Or, can use \nocite{wickham2016ggplot2, garnier2024viridis}. -->

The views expressed in this manuscript are those of the authors and do not necessarily represent the views or policies of USEPA, NOAA, or USFS. Any mention of trade names, products, or services does not imply an endorsement by the U.S. government, USEPA, NOAA, or USFS. USEPA, NOAA, or USFS do not endorse any commercial products, services or enterprises.

# References
