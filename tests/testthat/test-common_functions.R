#import data for testing
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



#tests for .bound function for bounding denominator of clever covariates
test_that(".bound works", {
  expect_equal(.bound(0.5,100), 0.5)
})

test_that("Bounds for denominator of clever covariate between 0 and 1", {
  expect_true(.bound(0.001,100)>0)
  expect_true(.bound(0.999,100)<1)
})


test_that("Error if denominator of clever covariate 0 for any observation", {
  expect_error(.bound(0,100), "Denominator of clever covariate 0 for at least one observation. Positivity violation?")
})

#tests for preprocess function
test_that("Message if removing observations missing treatment variable", {
  dat1 <- data
  dat1$treatment[1] <- NA
  expect_message(preprocess(txinrwd=TRUE, data=dat1, study="study", covariates=c("aged", "sex", "momedu", "hfiacat"), treatment_var="treatment", treatment=1, outcome="laz"), "Removing observations with missing treatment variable.")
})

test_that("Confirm no covariates in RWD not represented in RCT if txinrwd=FALSE (avoid positivity violation)", {
  dat1 <- data
  dat1 <- dat1[-which(dat1$study==2 & dat1$treatment==1),]
  check <- preprocess(txinrwd=FALSE, data=dat1, study="study", covariates=c("aged"), treatment_var="treatment", treatment=1, outcome="laz")
  expect_true(min(check$aged[which(check$S>1)]) >= min(check$aged[which(check$S==1)]))
  expect_true(max(check$aged[which(check$S>1)]) <= max(check$aged[which(check$S==1)]))
})

test_that("Confirm does not trim observations with NCO values outside of RCT if adjustnco==FALSE", {
  dat1 <- data
  dat1 <- dat1[-which(dat1$study==2 & dat1$treatment==1),]
  dat1[which(dat1$study==2),]$Nlt18scale[1] <- -2
  dat1[which(dat1$study==2),]$Nlt18scale[2] <- 10
  check <- preprocess(txinrwd=FALSE, data=dat1, study="study", covariates=c("aged"), treatment_var="treatment", treatment=1, outcome="laz", NCO="Nlt18scale", adjustnco = FALSE)
  expect_equal(min(check$nco[which(check$S>1)]), -2)
  expect_equal(max(check$nco[which(check$S>1)]), 10)
})

#tests for apply_selectorfunc function

#tests for validpreds function

#tests for limitdistvar function


