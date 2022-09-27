data(wash)
dat <- wash[which(wash$study %in% c(1,2)),]

library(SuperLearner)
txinrwd=TRUE
data=dat
study="study"
covariates=c("aged", "sex", "momedu", "hfiacat")
treatment_var="treatment"
treatment=1
outcome="laz"
NCO="Nlt18scale"
Delta=NULL
Delta_NCO=NULL
pRCT=0.5
V=10
Q.SL.library=c("SL.glm")
g.SL.library=c("SL.glm")
Q.discreteSL=TRUE
g.discreteSL=TRUE
family="gaussian"
family_nco="gaussian"
fluctuation = "logistic"
comparisons = list(c(1),c(1,2))
adjustnco = FALSE
target.gwt = TRUE

test_that("Y scaled appropriately if continuous outcome with logistic fluctuation", {

  out <- selector_func_notxrwd(train_s, data, Q.SL.library=c("SL.glm"), d.SL.library=NULL, g.SL.library=c("SL.glm"),
                               pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                               NCO=NULL, Delta=NULL, Delta_NCO = NULL,
                               adjustnco=FALSE, target.gwt=TRUE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE)
  expect_true(all(out$Y>=0))
  expect_true(all(out$Y<=1))
})
