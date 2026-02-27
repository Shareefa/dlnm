# Analysis Quirks and Discoveries

## mixmeta qtest() I-squared for Multivariate Models

The `mixmeta::qtest()` function does NOT return an I-squared field for multivariate models (multi-coefficient meta-analysis). I² must be computed manually:

```r
qt <- qtest(model)
# For multivariate models, qt$Q and qt$df are per-outcome
I2 <- pmax(0, (qt$Q - qt$df) / qt$Q * 100)
```

This was discovered during the setup-mixmeta-pipeline feature. The per-outcome Q statistics require manual I² computation as `max(0, (Q-df)/Q*100)`.

## crossbasis() NA Behavior for Incomplete Lag Windows

`crossbasis()` returns NA for the first `lag[2]` rows where the lag window is incomplete. For lag 0-21, the first 21 rows of each city/group will have NA values in the cross-basis matrix.

When generating synthetic truth data, set these NA effects to 0 (baseline) since they won't affect model fitting — the GLM will exclude NA rows from estimation.

```r
cb_effect <- as.numeric(unclass(cb) %*% beta_true)
cb_effect[is.na(cb_effect)] <- 0  # baseline for incomplete lag window
```

## Single-Stage vs Two-Stage CI Coverage

Single-stage DLNM (all cities in one GLM with group boundaries) produces **very low CI coverage** (~2%) because it computes conditional SEs that do not account for between-city heterogeneity. Two-stage DLNM+mixmeta produces much higher coverage (~76%) by properly modeling between-city variability through random effects.

This is **expected statistical behavior**, not a bug. The single-stage approach underestimates uncertainty because it treats city-specific effects as fixed, while the two-stage approach correctly pools heterogeneous estimates.
