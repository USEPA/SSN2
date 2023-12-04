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
  - name: Erin Peterson
    affiliation: 2
  - name: Jay M. Ver Hoef
    affiliation: 3
  - name: Alan Pearse
    affiliation: 4
  - name: Dan Isaak
    affiliation: 5
affiliations:
 - name: Office of Reserach and Development, United States Environmental Protection Agency
   index: 1
 - name: School of Mathematical Sciences, Queensland University of Technology, Brisbane, QLD Australia 4000
   index: 2
 - name: NMFS Alaska Fisheries Science Center, United States National Oceanic and Atmospheric Administration
   index: 3
 - name: NIASRA, School of Mathematics and Applied Statistics, University of Wollongong
   index: 4
 - name: Rocy Mountain Research Station, United States Forest Service
   index: 5
citation_author: Dumelle et. al.
date: 11 November 2023
year: 2023
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

The `SSN2` **R** package provides tools for spatial statistical modeling and prediction on stream (river) networks. `SSN2` is the successor to the `SSN` **R** package [@ver2014ssn], which was archived alongside broader changes in the **R**-spatial ecosystem [@nowosad2023] that included 1) the retirement of `rgdal` [@bivand2021rgdal], `rgeos` [@bivand2020rgeos], and `maptools` [@bivand2021maptools] and 2) the lack of active development of `sp` [@bivand2013applied]. `SSN2` leverages modern **R**-spatial tools like `sf` [@pebesma2018sf] and provides many useful modeling features that were not feasible to implement in `SSN`. 

# Statement of Need

Streams provide vital aquatic services that sustain wildlife, provide drinking and irrigation water, and support recreational and cultural activities.  Data are often collected at various locations on a stream network and used to characterize some scientific phenomenon in the stream. For example, a manger may need to know how the amount of a hazardous chemical changes throughout a stream network to inform mitigation efforts. Comprehensive formulations of SSN models are provided by @ver2010moving, @peterson2010mixed, and @ver2014ssn. 

SSN models use a spatial statistical modeling framework [@cressie1993statistics] to describe unique and complex dependencies on a stream network resulting from a branching network structure, directional water flow, and differences in flow volume. SSN models relate a response variable to one or more explanatory variables, a spatially independent error term (i.e., nugget), and up to three spatially dependent error terms: tail-down errors, tail-up errors, and Euclidean errors. Tail-down errors restrict spatial dependence to flow-connected sites (i.e., water flows from an upstream to a downstream site) and incorporate spatial weights (i.e., additive function) to describe the branching network between them. Tail-up errors describe spatial dependence between both flow-connected and flow-unconnected (i.e., sites that share a common downstream junction but not flow) sites, but spatial weights are not required. Euclidean errors describe spatial dependence between sites based on Euclidean distance and are governed by factors not confined to the stream network like regional geology. The length-scales of spatial dependence in the tail-up, tail-down, and Euclidean errors are controlled by separate range parameters. Next we show how to use the `SSN2` **R** package to fit and inspect SSN models and make predictions at unobserved locations on a stream network.

# Package Overview

Before fitting SSN models using `SSN2`, stream network data must be pre-processed using the STARS toolset for ArcGIS Desktop versions 9.3x - 10.8x [@peterson2014stars]. STARS is used to create a `.ssn` object (i.e., folder), which contains all the spatial, topological, and attribute information needed to fit stream network using `SSN2`. Shapefiles and text files residing in the `.ssn` object are read into **R** (using `ssn_import()`) and placed into a special list we call an SSN object. The SSN object contains geometry information, observed data, and data requiring prediction.

The `SSN2` packages comes with an example `.ssn` object that represents a stream network for the Middle Fork Basin of the Salmon River in Idaho, USA during 2004. To use this example data, we must copy the `.ssn` folder that comes shipped with `SSN2` to a temporary directory and store the temporary directory's file path:

```r
library(SSN2)
copy_lsn_to_temp()
path <- paste0(tempdir(), "/MiddleFork04.ssn")
```

Copying to the temporary directory is only used in this illustrative example, as users will read in the `.ssn` object from an appropriate file path stored on their computer (though file paths for SSN objects may be updated using `ssn_update_path()`). Next we import the observed data and prediction sites (`pred1km`):

```r
mf04p <- ssn_import(path, predpts = "pred1km")
```

We visualize the stream network, observed sites, and prediction sites using `ggplot2` [@wickham2016ggplot2] by running

```r
library(ggplot2)
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$preds$pred1km) +
  geom_sf(data = mf04p$obs, color = "brown", size = 1.5) +
  theme_bw()
```

\begin{figure}

{\centering \includegraphics[width=0.75\linewidth]{paper_files/figure-latex/steamnetwork-1} 

}

\caption{Middle Fork 2004 stream newtork. Observed sites are represnted by brown, closed circles. Prediction sites are represented by black, closed circles.}\label{fig:steamnetwork}
\end{figure}

We supplement the `.ssn` object with flow-connected and flow-unconnected distance matrices that are required for statistical modeling by running

```r
ssn_create_distmat(mf04p, predpts = "pred1km", overwrite = TRUE)
```

Suppose we model summer water temperature (`Summer_mn`) as a function of elevation (`ELEV_DEM`) and precipitation (`AREAWTMAP`) with a exponential, spherical, and Gaussian structures for the tail-up, tail-down, and Euclidean errors, respectively. We fit and summarize this model using `SSN2` by running

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

A summary of the fitted model looks similar to a summary returned by `lm()` but also returns spatial dependence parameter estimates:

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

`SSN2` leverages the `tidy()`, `glance()`, and `augment()` functions [@robinson2021broom] to tidy, glance at, and augment (with diagnostics) the fitted model:

```r
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

```r
glance(ssn_mod)
```

```
## # A tibble: 1 x 9
##       n     p  npar value   AIC  AICc logLik deviance pseudo.r.squared
##   <int> <dbl> <int> <dbl> <dbl> <dbl>  <dbl>    <dbl>            <dbl>
## 1    45     3     7  59.3  73.3  76.3  -29.6     41.9            0.612
```

```r
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

Specific helper functions (e.g., `coef()`, `AIC()`, `residuals()`) can be used to obtain the same quantities returned by `tidy()`, `glance()`, and `augment()`:

```r
coef(ssn_mod)
```

```
## (Intercept)    ELEV_DEM   AREAWTMAP 
## 76.19504087 -0.02690478 -0.00909941
```

```r
AIC(ssn_mod)
```

```
## [1] 73.2623
```

```r
head(residuals(ssn_mod))
```

```
##         1         2         3         4         5         6 
## -3.066413 -2.204147 -2.252004 -2.175337 -2.131527 -2.162417
```

Prediction at the prediction sites is performed using `augment` (or `predict()`):

```r
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
## 10    14.9   14.5   15.2 (-1529689 2537313)
## # i 165 more rows
```

Generalized linear models for binary, count, proportion, and skewed data are available via the `ssn_glm()` function. Simulating data on a stream network is performed via `ssn_simulate()`.

# Discussion

SSN models are invaluable tools for statistical analysis of stream network data and help to maintain and improve vital services that stream ecosystems provide. They have been employed to better understand an manage water quality [@scown2017improving; @mcmanus2020variation], ecosystem metabolism [@rodriguez2019estimating], and climate change impacts on freshwater ecosystems [@ruesch2012projected; @isaak2017norwest], as well as generate aquatic population estimates [@isaak2017scalable], inform conservation planning [@rodriguez2019spatial; @sharma2021dendritic], and assess restoration activities [@fuller2022riparian], among other applications.

There are several spatial modeling packages in **R**, including `geoR` [@ribiero2022geoR], `gstat` [@pebesma2004gstat], `FRK` [@sainsbury2002frk], `fields` [@nychka2021fields], and `R-INLA` [@lindgren2015bayesian], `spmodel` [@dumelle2023spmodel], among others. However, these packages fail to account for the intricacies of stream networks. `rtop`  [@skoien2014rtop] allows for spatial prediction on stream networks but fails to provide options for model fitting and diagnostics. Thus, `SSN2` is the most complete tool available in **R** for working with SSN models. To learn more about `SSN2`, visit our CRAN webpage at [https://CRAN.R-project.org/package=SSN2](https://CRAN.R-project.org/package=SSN2).

# Acknowledgements

The views expressed in this manuscript are those of the authors and do not necessarily represent the views or policies of USEPA, NOAA, or USFS. Any mention of trade names, products, or services does not imply an endorsement by the U.S. government, USEPA, NOAA, or USFS. USEPA, NOAA, or USFS do not endorse any commercial products, services or enterprises.

# References
