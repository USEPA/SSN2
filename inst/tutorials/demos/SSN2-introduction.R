#:::::::::::::::::::::::::::::::::::::
# Data Processing ####
#:::::::::::::::::::::::::::::::::::::


# Load packages ####
library(SSN2)
library(ggplot2)
library(sf)

# Copy Data Folder ####
## copy MiddleFork04 to the temporary directory
copy_lsn_to_temp()
## Find the path to the MiddleFork04 in the temporary directory
path <- paste0(tempdir(), "/MiddleFork04.ssn")

# Import MiddleFork04 SSN ####
## Import data and prediction points
mf04p <- ssn_import(
  path = path,
  predpts = c("pred1km", "CapeHorn"),
  overwrite = TRUE
)
## print a summary
summary(mf04p)
## print names of prediction data
names(mf04p$preds)
## visualize stream network, prediction data, observed data
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$preds$pred1km, pch = 17, color = "blue") +
  geom_sf(data = mf04p$obs, color = "brown", size = 2) +
  theme_bw()
## create distance matrices needed for spatial statistical modeling
ssn_create_distmat(
  ssn.object = mf04p,
  predpts = c("pred1km", "CapeHorn"),
  among_predpts = TRUE,
  overwrite = TRUE
)

#:::::::::::::::::::::::::::::::::::::
# Statistical Modeling ####
#:::::::::::::::::::::::::::::::::::::

# Exploratory Visualizations ####
## Visualize temperature spatially
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$obs, aes(color = Summer_mn), size = 2) +
  scale_color_viridis_c(limits = c(-1.5, 17), option = "H") +
  theme_bw()

## Visualize relationship between elevation and temperature
ggplot(data = mf04p$obs, aes(x = ELEV_DEM, y = Summer_mn)) +
  geom_point(size = 2) +
  theme_bw()

## Visualize relationship between precipitation and temperature
ggplot(data = mf04p$obs, aes(x = AREAWTMAP, y = Summer_mn)) +
  geom_point(size = 2) +
  theme_bw()

## Create Torgegram
tg <- Torgegram(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  type = c("flowcon", "flowuncon", "euclid")
)
## Visualize Torgegram
plot(tg)

# Model Building ####
## Create a model with all possible covariance components
ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)
## Summarize the model
summary(ssn_mod)

## Inspect the variance components
varcomp(ssn_mod)

## Glance at the model fit
glance(ssn_mod)

## Fit a model with just taildown and nugget (notice no additive needed)
ssn_mod2 <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  taildown_type = "spherical"
)

## Fit a model with just nugget (equivalent to lm)
ssn_mod3 <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p
)

## Glance at all model fits (look for lowest AIC)
glances(ssn_mod, ssn_mod2, ssn_mod3)

## leave-one-out cross validation for each model; find RMSPE
loocv_mod <- loocv(ssn_mod)
loocv_mod
loocv_mod$RMSPE
loocv_mod2 <- loocv(ssn_mod2)
loocv_mod2$RMSPE
loocv_mod3 <- loocv(ssn_mod3)
loocv_mod3$RMSPE

## Fit models with maximum likelihood (this allows us to compare models
## with different fixed effect structures and covariance structures)
ml_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea",
  estmethod = "ml"
)
## Model without ELEV_DEM
ml_mod2 <- ssn_lm(
  formula = Summer_mn ~ AREAWTMAP,
  ssn.object = mf04p,
  euclid_type = "gaussian",
  estmethod = "ml"
)

## glance at model fit and leave-one-out
glances(ml_mod, ml_mod2)
loocv_mod_ml <- loocv(ml_mod)
loocv_mod_ml$RMSPE
loocv_mod_ml2 <- loocv(ml_mod2)
loocv_mod_ml2$RMSPE

## Refit final model using REML
loocv_mod$RMSPE
loocv_mod_ml$RMSPE


#:::::::::::::::::::::::::::::::::::::
# Model Inspection ####
#:::::::::::::::::::::::::::::::::::::

## Confidence intervals on parameters
tidy(ssn_mod, conf.int = TRUE)

## Augment the data with diagnostics
aug_ssn_mod <- augment(ssn_mod)
aug_ssn_mod

## Write out the augmented data to a geopackage
st_write(aug_ssn_mod, paste0(tempdir(), "/aug_ssn_mod.gpkg"))

## Model visualization (fitted vs residuals [standardized])
plot(ssn_mod, which = 1)


#:::::::::::::::::::::::::::::::::::::
# Prediction ####
#:::::::::::::::::::::::::::::::::::::

## Make predictions at all locations in the pred1km data
predict(ssn_mod, newdata = "pred1km")

## Augment prediction data with predictions
aug_preds <- augment(ssn_mod, newdata = "pred1km")
aug_preds[, ".fitted"]

## Visualize predictions ont he network
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = aug_preds, aes(color = .fitted), size = 2) +
  scale_color_viridis_c(limits = c(-1.5, 17), option = "H") +
  theme_bw()

## Write out augmented prediction data to a geopackage
st_write(aug_preds, paste0(tempdir(), "/aug_preds.gpkg"))

## Predict for all prediction data
predict(ssn_mod)
predict(ssn_mod, newdata = "all")

## Block prediction for the average in an area (here, CapeHorn)
predict(ssn_mod, newdata = "CapeHorn", block = TRUE, interval = "prediction")

#:::::::::::::::::::::::::::::::::::::
# Advanced Features ####
#:::::::::::::::::::::::::::::::::::::

# Specify initial values ####

## Fix de at one and fit model
euclid_init <- euclid_initial("gaussian", de = 1, known = "de")
euclid_init

ssn_init <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_initial = euclid_init,
  additive = "afvArea"
)
ssn_init

# Random Effects ####

## Random effect for network ID (two networks)
ssn_rand <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea",
  random = ~ as.factor(netID)
)
ssn_rand

# Partition Factors ####

## Partition Factor for network ID (even Euclidean component is independent
## across networks when specifying this)
ssn_part <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea",
  partition_factor = ~ as.factor(netID)
)
ssn_part

# Generalized Linear Models

## Visualize C16, the number of days the temperature exceeded 16 C
## This is a count variable
ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$obs, aes(color = C16), size = 2) +
  scale_color_viridis_c(option = "H") +
  theme_bw()

## Fit a poisson regression model with tailup and taildown components
ssn_pois <- ssn_glm(
  formula = C16 ~ ELEV_DEM + AREAWTMAP,
  family = "poisson",
  ssn.object = mf04p,
  tailup_type = "epa",
  taildown_type = "mariah",
  additive = "afvArea"
)

## Sumamrize the model fit
summary(ssn_pois)

## Fit a negative binomial model with tailup and taildown components
ssn_nb <- ssn_glm(
  formula = C16 ~ ELEV_DEM + AREAWTMAP,
  family = "nbinomial",
  ssn.object = mf04p,
  tailup_type = "epa",
  taildown_type = "mariah",
  additive = "afvArea"
)

# Compare via AIC, loocv
glances(ssn_pois, ssn_nb)
loocv_pois <- loocv(ssn_pois)
loocv_pois$RMSPE
loocv_nb <- loocv(ssn_nb)
loocv_nb$RMSPE

# Simulating Data on a Stream Network

## Specify covariance parameters
tu_params <- tailup_params("exponential", de = 0.4, range = 1e5)
td_params <- taildown_params("spherical", de = 0.1, range = 1e6)
euc_params <- euclid_params("gaussian", de = 0.2, range = 1e3)
nug_params <- nugget_params("nugget", nugget = 0.1)

## Set a reproducible seed
set.seed(2)

## Simulate Gaussian data
sims <- ssn_simulate(
  family = "gaussian",
  ssn.object = mf04p,
  network = "obs",
  additive = "afvArea",
  tailup_params = tu_params,
  taildown_params = td_params,
  euclid_params = euc_params,
  nugget_params = nug_params,
  mean = 0,
  samples = 1
)
## Inspect first few simulations
head(sims)

## Simulate binomial data
sims <- ssn_simulate(
  family = "binomial",
  ssn.object = mf04p,
  network = "obs",
  additive = "afvArea",
  tailup_params = tu_params,
  taildown_params = td_params,
  euclid_params = euc_params,
  nugget_params = nug_params,
  mean = 0,
  samples = 2
)
## Inspect first few simulations
head(sims)
