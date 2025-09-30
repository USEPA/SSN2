# SSN2 0.3.1

## Bug fixes

* Added a dependency on `spmodel` version at least `0.7.0` so that `AUROC()` methods are properly called.
* Fixed a bug introduced in version 0.3.0 that prevented `partition_factor` from working properly when non-`NULL` ([#32](https://github.com/USEPA/SSN2/issues/32)).

# SSN2 0.3.0

## Major Updates

* Added support for the `emmeans` **R** package for estimating marginal means of `ssn_lm()` and `ssn_glm()` models.
* Added support for applications to large data sets. The `ssn_create_bigdist()` function was added to create large distance matrices using the `filematrix` **R** package. Estimation for large data sets is performed by leveraging the `local` argument to `ssn_lm()` and `ssn_glm()`. Prediction for large data sets is performed by leveraging the `local` argument to `predict()` (and `augment()`). When `local` is used, `SSN2` looks for distance matrices created using `ssn_create_bigdist()`.
* Added support for Gaussian tail-up and tail-down covariance functions.


## Minor Updates

* Updated `ssn_import()` so that it does not force an overwrite of the `netgeom` column when it already exists.
* Add a `verbose` argument to `ssn_import()`, `ssn_import_predpts()`, and `createBinaryID()` to control whether warning messages are printed to the **R** console.
* Added the `na.action` argument to `predict.ssn_lm()` and `predict.ssn_glm()` functions to clarify that missing values in `newdata` return an error.
* Changed the `type` argument in `augment()` for `ssn_glm()` models to `type.predict` to match `broom::augment.glm()`.
* `augment()` for `ssn_glm()` models now returns fitted values on the link scale by default to match `broom::augment.glm()`.
* Added a `type.residuals` argument for `ssn_glm()` models to match `broom::augment.glm()`.
* Updated `logLik()` to match `lm()` and `glm()` behavior. `logLik()` now returns a vector with class `logLik` and attributes `nobs` and `df`.
* Added support for using `AIC()` and `BIC()` from `stats` and removed `SSN2`-specific `AIC()` methods.
* Added a `warning` argument to `glances()` that determines whether relevant warnings should be displayed or not.
* Added a warning message to `glances()` about interpreting likelihood-based statistics (e.g., AIC, AICc, BIC) when a one model has `estmethod = "ml"` and another model has `estmethod = "reml"`.
* Added a warning message to `glances()` about interpreting likelihood-based statistics (e.g., AIC, AICc, BIC) when two models with `estmethod = "reml"` have distinct `formula` arguments.
* Added a warning message to `glances()` about interpreting likelihood-based statistics (e.g., AIC, AICc, BIC) when two models have different sample sizes.
* Added a warning message to `glances()` about interpreting likelihood-based statistics (e.g., AIC, AICc, BIC) when two models have different family supports (which can happen with `ssn_glm()` models).
* Added a `cloud` argument to `Torgegram()` to return a cloud Torgegram.
* Added the ability to pass custom `cex` to `plot.Torgegram()`.
* Added a robust semivariogram option to `Torgegram()`; see the `robust`argument to `Torgegram()`
* Added an `AUROC()`  function to compute the area under the receiver operating characteristic (AUROC) curve for `ssn_glm` models when `family` is `"binomial"` and the response is binary (i.e., represents a single success or failure).
* Added a `type` argument to `loocv()` when `cv_predict = TRUE` and using `ssn_glm()` models so that predictions may be obtained on the link or response scale.
* Added support for `"terms"` prediction for `ssn_lm()` and `ssn_glm()` models.
* Added `scale` and `df` arguments to `predict()` for `ssn_lm()` models.
* Add `dispersion` argument to `predict()` for `ssn_glm()` models.
* Minor documentation updates.
* Minor code maintenance updates.

## Bug Fixes

* Fixed a bug that caused incorrect degrees of freedom for the likelihood ratio test (`anova(model1, model2)`) when `estmethod` is `"ml"` for both models [#25](https://github.com/USEPA/SSN2/issues/25).
* Fixed a bug that caused an error in `anova(object1, object2)` when the name of `object1` had special characters (e.g., `$`).

# SSN2 0.2.1

## Minor Updates

* Enhanced numeric stability of deviance and pseudo R-squared for `ssn_glm()` models when `family = "beta"` [(#23)](https://github.com/USEPA/SSN2/issues/23).
* Updated `reexport.Rd` to reflect changes in `spmodel v0.8.0`'s handling of `AIC()` and `AICc()`.

# SSN2 0.2.0

## Major Updates

* Significant testing, documentation, and auxiliary (e.g., `README.md`) updates as part of a submission to *Journal of Open Source Software*. Relevant issues associated with the review are available at [#11](https://github.com/USEPA/SSN2/issues/11), [#12](https://github.com/USEPA/SSN2/issues/12), [#13](https://github.com/USEPA/SSN2/issues/13), [#14](https://github.com/USEPA/SSN2/issues/14), [#15](https://github.com/USEPA/SSN2/issues/15), [#16](https://github.com/USEPA/SSN2/issues/16), [#17](https://github.com/USEPA/SSN2/issues/17), [#20](https://github.com/USEPA/SSN2/issues/20), [#21](https://github.com/USEPA/SSN2/issues/21). The review is [linked here](https://github.com/openjournals/joss-reviews/issues/6389).
* Added support for geopackage file formats in the `.ssn` folder that is accessed when importing SSN objects via `ssn_import()`.

## Minor Updates

* Added `ssn_names()` to return column names in the `edges`, `obs`, and `preds` elements of an SSN object.
* Changed `Matrix::rankMatrix(X, method = "tolNorm2")` to `Matrix::rankMatrix(X, method = "qr")` to enhance stability when determining linear independence in `X`, the design matrix of explanatory variables.
* Replaced an error message with a warning message when `X` has perfect collinearities (i.e., is not full rank).
* Removed `format_additive` argument from `ssn_import()` because of transition to geopackage support, which eliminates the need to convert additive function values to text.
* Added the `create_netgeom()` function to create the network geometry column for the `edges`, `obs`, and `preds` elements in an SSN object.
* Minor vignette updates.
* Minor documentation updates.

## Bug Fixes

* Fixed a bug in `SSN_to_SSN2()` that caused an error using `ssn_write()` with no prediction sites.
* Replaced `names.SSN()` with `ssn_names()`, as `names.SSN()` prevented proper naming of elements in the SSN object.

# SSN2 0.1.1

## Minor Updates

* Changed network geometry name from `netgeometry` to `netgeom` to avoid exceeding the 10 character limit for column/field names while writing to shapefiles [(#2)](https://github.com/USEPA/SSN2/issues/2).
* Added an error message when `family` is missing in `ssn_glm()` [(#8)](https://github.com/USEPA/SSN2/issues/8).
* Added a deprecation warning for `SSN_to_SSN2()`.
* Minor stability updates.
* Minor error message updates.
* Minor documentation updates.

## Bug Fixes

* Fixed a bug in `Torgegram()` that prevented intended computation when `cutoff` was specified.
* Fixed a bug in `plot.Torgegram()` that occasionally prevented proper spacing of the legend.
* Fixed a bug that prevented proper printing of the dispersion parameter from `ssn_glm()` model objects (and their summaries) when all covariance parameters were known.
* Fixed a bug that prevented simulation when `euclid_type` was `"none"`.
* Fixed a bug that could cause improper prediction behavior when `taildown_type` was `"spherical"`.
* Fixed a bug that printed response residuals instead of deviance residuals for `ssn_glm()` objects.

# SSN2 0.1.0

* Initial CRAN submission.
