# SSN2 0.1.2

## Minor Updates

* Added `ssn_names()` to return column names in the `edges`, `obs`, and `preds` elements of an SSN object.

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
