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
