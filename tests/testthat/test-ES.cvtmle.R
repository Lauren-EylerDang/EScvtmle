testthat::test_that("errors", {
  testthat::expect_error(
    ES.cvtmle(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, pRCT, V=10, Q.SL.library, d.SL.library, g.SL.library, Q.discreteSL, d.discreteSL, g.discreteSL, family, family_nco, fluctuation = "logistic", comparisons = list(c(1),c(1,2),c(1,2,3)), adjustnco = FALSE, target.gwt = TRUE),
    "Package currently compares two experiments. Check back for updates to compare multiple experiments."
  )
})
