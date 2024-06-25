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
 - name: EP Consulting and Queensland University of Technology, Centre for Data Science Brisbane Australia 4000 
   index: 2
 - name: NMFS Alaska Fisheries Science Center, United States National Oceanic and Atmospheric Administration (NOAA)
   index: 3
 - name: NIASRA, School of Mathematics and Applied Statistics, University of Wollongong
   index: 4
 - name: Rocky Mountain Research Station, United States Forest Service (USFS)
   index: 5
citation_author: Dumelle et. al.
date: 11 November 2023
year: 2023
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
editor_options: 
  chunk_output_type: console
---

# Summary

The `SSN2` **R** package provides tools for spatial statistical modeling, parameter estimation, and prediction on stream (river) networks. `SSN2` is the successor to the `SSN` **R** package [@ver2014ssn], which was archived alongside broader changes in the **R**-spatial ecosystem [@nowosad2023] that included 1) the retirement of `rgdal` [@bivand2021rgdal], `rgeos` [@bivand2020rgeos], and `maptools` [@bivand2021maptools] and 2) the lack of active development of `sp` [@bivand2013applied]. `SSN2` maintains compatibility with the input data file structures used by `SSN` but leverages modern **R**-spatial tools like `sf` [@pebesma2018sf] and provides many useful features that were not available in `SSN`, including new modeling and helper functions, updated fitting algorithms, and simplified syntax consistent with other **R** generic functions.

# Statement of Need

Streams provide vital aquatic services that sustain wildlife, provide drinking and irrigation water, and support recreational and cultural activities.  Data are often collected at various locations on a stream network and used to characterize spatial patterns in stream phenomena. For example, a manager may need to know how the amount of a hazardous chemical changes throughout a stream network to inform mitigation efforts. Comprehensive formulations of spatial stream network (SSN) models are provided by @ver2010moving, @peterson2010mixed, and @ver2014ssn. The `SSN2` **R** package is designed to help users fit SSN models to their stream network data.

SSN models use a spatial statistical modeling framework [@cressie1993statistics] to describe unique and complex dependencies on a stream network resulting from a branching network structure, directional water flow, and differences in flow volume. SSN models relate a continuous or discrete response variable to one or more explanatory variables, a spatially independent error term (i.e., nugget), and up to three spatially dependent error terms: tail-up errors, tail-down errors, and Euclidean errors. Tail-up errors restrict spatial dependence to flow-connected sites (i.e., water flows from an upstream to a downstream site) and incorporate spatial weights through an additive function to describe the branching network between sites. Tail-down errors describe spatial dependence between both flow-connected and flow-unconnected (i.e., sites that share a common downstream junction but not flow) sites, but spatial weights are not required. Euclidean errors describe spatial dependence between sites based on straight-line distance and are governed by factors not confined to the stream network like regional geology. The length-scales of spatial dependence in the tail-up, tail-down, and Euclidean errors are controlled by separate range parameters. In this paper, we show how to use the `SSN2` **R** package to fit SSN models, inspect SSN models, and use SSN models to make predictions at unobserved locations on a stream network. 

# Package Overview

Before fitting SSN models using `SSN2`, stream network and observation data sets must be pre-processed either by using the STARS toolset for ArcGIS Desktop versions 9.3x - 10.8x [@peterson2014stars], by using the `openSTARS` **R** package [@kattwinkel2020preparing], which leverages open-source GRASS GIS, or by using the `SSNbler` **R** package [@peterson2024SSNbler], a new, R-based version of STARS that is available on GitHub, will soon be available on CRAN, and contains several useful resources that guide users through these pre-processing steps. Pre-processing using either STARS, `openSTARS`, or `SSNbler` ends with the creation of a `.ssn` folder, which is non-proprietary and has all the spatial, topological, and attribute information needed to fit models to data on a stream network using `SSN2`. Relevant files residing in the `.ssn` folder are read into **R** (using `ssn_import()`) and placed into a list and called an SSN object. The SSN object contains geometry and topological information about the stream reaches and sites, as well as observed data and data for prediction at unsampled sites. 

`SSN2` is first installed from CRAN:

```r
install.packages("SSN2")
```

Then, `SSN2` is loaded into our current **R** session:

```r
library(SSN2)
```

The `SSN2` packages comes with an example `.ssn` folder called `MiddleFork04.ssn` that represents water temperatures recorded from a stream network in the Middle Fork of the Salmon River in Idaho, USA during 2004. 

Several functions in `SSN2` for reading and writing data (which we use shortly) directly manipulate the `.ssn` folder. As to avoid directly manipulating the `MiddleFork04.ssn` data installed alongside `SSN2`, `MiddleFork04.ssn` is instead be copied it into a temporary directory and the relevant path to directory stored:

```r
copy_lsn_to_temp()
path <- paste0(tempdir(), "/MiddleFork04.ssn")
```

The `copy_lsn_to_temp()` function is only used when working with `MiddleFork04.ssn` and generally, `path` should indicate a permanent directory on your machine that points towards your `.ssn` object. After specifying `path`, the stream reaches, observed sites, and prediction sites (`pred1km`) are imported:

```r
mf04p <- ssn_import(path, predpts = "pred1km")
```

The stream network, observed sites, and prediction sites (Figure$~$\ref{fig:steam-network})  are visualized [@wickham2016ggplot2]:

```r
library(ggplot2)
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$preds$pred1km, pch = 17, color = "blue") +
  geom_sf(data = mf04p$obs, color = "brown", size = 2) +
  theme_bw()
```

\begin{figure}

{\centering \includegraphics[width=0.75\linewidth]{paper_files/figure-latex/steam-network-1} 

}

\caption{Middle Fork 2004 stream networks. Observed sites are represented by brown, closed circles. Prediction sites are represented by blue, closed triangles.}\label{fig:steam-network}
\end{figure}

Prior to statistical modeling, the `.ssn` object must be supplemented with hydrologic distance matrices [@ver2010moving]:

```r
ssn_create_distmat(mf04p, predpts = "pred1km", overwrite = TRUE)
```

Of particular interest in this example is summer mean stream temperature (`Summer_mn`) in degrees Celsius, visualized [@garnier2024viridis] via :

```r
ggplot() +
    geom_sf(data = mf04p$edges) +
    geom_sf(data = mf04p$obs, aes(color = Summer_mn), size = 2) +
    scale_color_viridis_c(name = "Obs.", option = "H", limits = c(-3.5, 17)) +
    theme_bw()
```

\begin{figure}

{\centering \includegraphics[width=0.75\linewidth]{paper_files/figure-latex/steam-obs-1} 

}

\caption{Observed Middle Fork 2004 mean summer temperatures (Celsius). }\label{fig:steam-obs}
\end{figure}


Wee will model summer mean stream temperature as a function of elevation (`ELEV_DEM`) and watershed-averaged precipitation (`AREAWTMAP`) with exponential, spherical, and Gaussian structures for the tail-up, tail-down, and Euclidean errors, respectively. Using `ssn_lm()`, the model is fit:

```r
ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)
```

The `ssn_lm()` function is designed to be similar in syntax and structure to the `lm()` function in base **R** for fitting nonspatial linear models. Additionally, `SSN2` accommodates various S3 methods for commonly-used **R** generic functions that operate on model objects. For example, the generic function `summary()` is used to summarize the fitted model:

```r
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
##      Min       1Q   Median       3Q      Max 
## -2.73430 -1.43161 -0.04368  0.83251  1.39377 
## 
## Coefficients (fixed):
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept) 78.214857  12.189379   6.417 1.39e-10 ***
## ELEV_DEM    -0.028758   0.005808  -4.952 7.35e-07 ***
## AREAWTMAP   -0.008067   0.004125  -1.955   0.0505 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Pseudo R-squared: 0.4157
## 
## Coefficients (covariance):
##               Effect     Parameter   Estimate
##   tailup exponential  de (parsill)  1.348e+00
##   tailup exponential         range  8.987e+05
##   taildown spherical  de (parsill)  2.647e+00
##   taildown spherical         range  1.960e+05
##      euclid gaussian  de (parsill)  1.092e-04
##      euclid gaussian         range  1.805e+05
##               nugget        nugget  1.660e-02
```

`SSN2` leverages the `tidy()`, `glance()`, and `augment()` generic functions [@robinson2021broom] to inspect the fitted model and provide diagnostics:

```r
tidy(ssn_mod, conf.int = TRUE)
```

```
## # A tibble: 3 x 7
##   term        estimate std.error statistic  p.value conf.low   conf.high
##   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>       <dbl>
## 1 (Intercept) 78.2      12.2          6.42 1.39e-10  54.3    102.       
## 2 AREAWTMAP   -0.00807   0.00413     -1.96 5.05e- 2  -0.0162   0.0000187
## 3 ELEV_DEM    -0.0288    0.00581     -4.95 7.35e- 7  -0.0401  -0.0174
```

```r
glance(ssn_mod)
```

```
## # A tibble: 1 x 9
##       n     p  npar value   AIC  AICc logLik deviance pseudo.r.squared
##   <int> <dbl> <int> <dbl> <dbl> <dbl>  <dbl>    <dbl>            <dbl>
## 1    45     3     7  76.6  90.6  93.7  -38.3     41.8            0.416
```

```r
aug_mod <- augment(ssn_mod)
subset(aug_mod, select = c(Summer_mn, .fitted, .resid, .hat, .cooksd))
```

```
## Simple feature collection with 45 features and 5 fields
## Geometry type: POINT
## Dimension:     XY
## Bounding box:  xmin: -1530805 ymin: 920324.3 xmax: -1503079 ymax: 931036.6
## Projected CRS: USA_Contiguous_Albers_Equal_Area_Conic
## # A tibble: 45 x 6
##    Summer_mn .fitted  .resid    .hat   .cooksd            geometry
##        <dbl>   <dbl>   <dbl>   <dbl>     <dbl>         <POINT [m]>
##  1      14.9    14.1  0.770  0.0724  0.00274   (-1528194 929550.4)
##  2      14.7    14.0  0.714  0.0569  0.0000449 (-1528222 928237.7)
##  3      14.6    13.8  0.776  0.0629  0.00259   (-1528485 927846.1)
##  4      15.2    14.8  0.427  0.125   0.0471    (-1519790 930112.1)
##  5      14.5    14.5 -0.0437 0.0359  0.0343      (-1520336 929772)
##  6      15.3    14.3  1.01   0.0220  0.00329   (-1524599 930808.7)
##  7      15.1    14.3  0.797  0.0178  0.000105  (-1525729 930933.4)
##  8      14.9    14.1  0.833  0.00213 0.0000813 (-1527966 929774.7)
##  9      15.0    13.9  1.06   0.0560  0.000182  (-1528257 929648.2)
## 10      15.0    13.9  1.15   0.0471  0.00684   (-1528428 929476.2)
## # i 35 more rows
```

Specific generic helper functions (e.g., `coef()`, `AIC()`, `residuals()`) can be used to obtain the same quantities returned by `tidy()`, `glance()`, and `augment()`:

```r
coef(ssn_mod)
```

```
##  (Intercept)     ELEV_DEM    AREAWTMAP 
## 78.214856580 -0.028758302 -0.008066962
```

```r
AIC(ssn_mod)
```

```
## [1] 90.63532
```

```r
head(residuals(ssn_mod))
```

```
##           1           2           3           4           5           6 
##  0.77010720  0.71389871  0.77644852  0.42749165 -0.04368363  1.00936419
```

Spatial prediction (i.e., Kriging) at the unobserved sites is performed using the generic functions `predict()` or `augment()`:

```r
aug_pred <- augment(ssn_mod, newdata = "pred1km", interval = "prediction")
subset(aug_pred, select = c(.fitted, .lower, .upper))
```

```
## Simple feature collection with 175 features and 3 fields
## Geometry type: POINT
## Dimension:     XY
## Bounding box:  xmin: -1530631 ymin: 914920.7 xmax: -1500020 ymax: 933466.4
## Projected CRS: USA_Contiguous_Albers_Equal_Area_Conic
## # A tibble: 175 x 4
##    .fitted .lower .upper            geometry
##      <dbl>  <dbl>  <dbl>         <POINT [m]>
##  1    14.7   14.4   15.0 (-1528406 928161.4)
##  2    14.7   14.3   15.2 (-1528202 928821.1)
##  3    14.9   14.5   15.2 (-1528173 929414.9)
##  4    14.4   12.6   16.1 (-1530218 926538.7)
##  5    14.5   12.7   16.2 (-1529466 926808.1)
##  6    14.5   14.1   14.9 (-1520657 929871.1)
##  7    15.0   14.7   15.4 (-1519866 930025.5)
##  8    14.7   13.7   15.6 (-1521823 930124.7)
##  9    14.9   13.8   16.0 (-1523183 930469.7)
## 10    15.2   14.7   15.8 (-1523860 930665.8)
## # i 165 more rows
```

Here, `.fitted` are the predictions, `.lower` are the lower bounds of 95% prediction intervals, and `.upper` are the upper bounds of 95% prediction intervals. Utilizing `augment()` makes the predictions straightforward to visualize:


```r
ggplot() +
    geom_sf(data = mf04p$edges) +
    geom_sf(data = aug_pred, aes(color = .fitted), size = 2) +
    scale_color_viridis_c(name = "Pred.", option = "H", limits = c(-3.5, 17)) +
    theme_bw()
```

\begin{figure}

{\centering \includegraphics[width=0.75\linewidth]{paper_files/figure-latex/steam-preds-1} 

}

\caption{Predicted Middle Fork 2004 mean summer temperatures (Celsius) spaced one kilometer apart. }\label{fig:steam-preds}
\end{figure}

Generalized spatial linear models for binary, count, proportion, and skewed data are available via the `ssn_glm()` function. Simulating data on a stream network is performed via `ssn_simulate()`.

# Discussion

SSN models are valuable tools for statistical analysis of data collected on stream networks and help improve inference about vital stream ecosystems. These models have been employed (using `SSN`) to better understand and manage water quality [@scown2017improving; @mcmanus2020variation], ecosystem metabolism [@rodriguez2019estimating], and climate change impacts on freshwater ecosystems [@ruesch2012projected; @isaak2017norwest], as well as generate aquatic population estimates [@isaak2017scalable], inform conservation planning [@rodriguez2019spatial; @sharma2021dendritic], and assess restoration activities [@fuller2022riparian], among other applications. The breadth and applicability of SSN models are further enhanced by data aggregation tools like the National Hydrography Dataset [@mckay2012nhdplus], National Stream Internet Project [@nagel2015national] and StreamCat [@hill2016stream].

There are several spatial modeling packages in **R**, including `geoR` [@ribiero2022geoR], `gstat` [@pebesma2004gstat], `FRK` [@sainsbury2002frk], `fields` [@nychka2021fields], `R-INLA` [@lindgren2015bayesian], and `spmodel` [@dumelle2023spmodel], among others. However, these packages fail to account for the intricacies of stream networks. `rtop`  [@skoien2014rtop] allows for spatial prediction on stream networks but fails to provide options for model fitting and diagnostics. Thus, `SSN2` is the most complete tool available in **R** for working with SSN models. To learn more about `SSN2`, visit our CRAN webpage at [https://CRAN.R-project.org/package=SSN2](https://CRAN.R-project.org/package=SSN2).

# Acknowledgements

The views expressed in this manuscript are those of the authors and do not necessarily represent the views or policies of USEPA, NOAA, or USFS. Any mention of trade names, products, or services does not imply an endorsement by the U.S. government, USEPA, NOAA, or USFS. USEPA, NOAA, or USFS do not endorse any commercial products, services or enterprises.

# References
