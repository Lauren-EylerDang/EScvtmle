#import data for testing
data(wash)
dat <- wash[which(wash$study %in% c(1,2)),]

library(SuperLearner)
txinrwd=TRUE
data=dat
study="study"
covariates=c("aged", "sex", "momedu")
treatment_var="treatment"
treatment=1
outcome="laz"
NCO="Nlt18scale"
Delta=NULL
Delta_NCO=NULL
pRCT=0.5
V=10
Q.SL.library=c("SL.glm")
d.SL.library=c("SL.glm")
g.SL.library=c("SL.glm")
Q.discreteSL=TRUE
d.discreteSL=TRUE
g.discreteSL=TRUE
family="gaussian"
family_nco="gaussian"
fluctuation = "logistic"
comparisons = list(c(1),c(1,2))
adjustnco = FALSE
target.gwt = TRUE
bounds = NULL



#tests for .bound function for bounding denominator of clever covariates
test_that(".bound works", {
  expect_equal(.bound(0.5,100), 0.5)
})

test_that("Bounds for denominator of clever covariate between 0 and 1", {
  expect_true(.bound(0.001,100)>0)
  expect_true(.bound(0.999,100)<1)
})

test_that(".bound respects user-specified bounds", {
  expect_equal(.bound(0.001,10000), 0.005428681)
  expect_true(.bound(0.001,10000, bounds = c(0.025,0.975))==0.025)
  expect_true(.bound(0.999,10000, bounds = c(0.025,0.975))==0.975)
})

#tests for preprocess function
test_that("Message if removing observations missing treatment variable", {
  dat1 <- data
  dat1$treatment[1] <- NA
  expect_message(preprocess(txinrwd=TRUE, data=dat1, study="study", covariates=c("aged", "sex", "momedu", "hfiacat"), treatment_var="treatment", treatment=1, outcome="laz"), "Removing observations with missing treatment variable.")
})

test_that("Confirm columns renamed appropriately", {
  check <- preprocess(txinrwd=TRUE, data=data, study="study", covariates=c("aged"), treatment_var="treatment", treatment=1, outcome="laz", NCO="Nlt18scale", adjustnco = FALSE)
  expect_true(all(c("S", "A", "Y", "nco") %in% colnames(check)))
})

test_that("Confirm treatment variable renamed appropriately", {
  dat1 <- data
  dat1[which(dat1$treatment==1),]$treatment <- "Treat"
  dat1[which(dat1$treatment==0),]$treatment <- "Ctrl"
  check <- preprocess(txinrwd=TRUE, data=dat1, study="study", covariates=c("aged"), treatment_var="treatment", treatment="Treat", outcome="laz", NCO="Nlt18scale", adjustnco = FALSE)
  expect_true(all(check$A %in% c(0,1)))
  expect_equal(length(which(dat1$treatment=="Treat")), length(which(check$A==1)))
})

test_that("Confirm no covariates in RWD not represented in RCT if txinrwd=FALSE (avoid positivity violation)", {
  dat1 <- data
  dat1 <- dat1[-which(dat1$study==2 & dat1$treatment==1),]
  dat1$momedu <- as.factor(dat1$momedu)
  dat1$momedu[which(dat1$momedu=="No education" & dat1$study==1)] <- "Primary (1-5y)"
  check <- preprocess(txinrwd=FALSE, data=dat1, study="study", covariates=c("aged"), treatment_var="treatment", treatment=1, outcome="laz")
  expect_true(min(check$aged[which(check$S>1)]) >= min(check$aged[which(check$S==1)]))
  expect_true(max(check$aged[which(check$S>1)]) <= max(check$aged[which(check$S==1)]))
  expect_true(length(which(check$momedu=="No education"))==0)
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

test_that("Confirm correct processing of missing outcomes for outcome and nco", {
  dat1 <- data
  dat1$laz[1:3] <- NA
  dat1$Nlt18scale[4:6] <- NA
  check <- preprocess(txinrwd=TRUE, data=dat1, study="study", covariates=c("aged", "sex"), treatment_var="treatment", treatment=1, outcome="laz", NCO="Nlt18scale", adjustnco = FALSE)
  expect_equal(check$NCO_delta[4:6], c(0,0,0))
  expect_equal(check$Delta[1:3], c(0,0,0))
})

#tests for apply_selectorfunc function
test_that("Confirm correct processing of missing outcomes for outcome and nco", {
  dat1 <- data
  check <- preprocess(txinrwd=TRUE, data=dat1, study="study", covariates=c("aged", "sex"), treatment_var="treatment", treatment=1, outcome="laz", NCO="Nlt18scale", adjustnco = FALSE)
  check2 <- apply_selector_func(txinrwd, train=check, data=check, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO=NULL, Delta=NULL, Delta_NCO=NULL, adjustnco=adjustnco, target.gwt=target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL, comparisons, bounds)
  expect_equal(length(check$Y[which(check$S==1)]), length(check2[[1]]$Y))
})


#tests for validpreds function
data$laz[sample(nrow(data),100)] <- NA
data <- preprocess(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO, Delta, Delta_NCO, adjustnco)
if("Delta" %in% colnames(data)){
  Delta = "Delta"
}

if("NCO_delta" %in% colnames(data)){
  Delta_NCO = "NCO_delta"
}
#Create cross-validation folds that preserve proportion of RCT (if txinwrd=TRUE) or of RCT controls (if txinrwd=FALSE) in validation sets
ids <- data$S

if(txinrwd==TRUE){
  ids[which(data$S==1)] <- 0
} else {
  ids[which(data$S==1 & data$A==0)] <- 0
}

folds <- make_folds(data, fold_fun = folds_vfold, V=V, strata_ids = ids)
data$v <- rep(NA, nrow(data))
for(v in 1:V){
  data$v[folds[[v]]$validation_set]<-v
}

results <- list()
selector <- list()
valid_initial <- list()

lambdatilde <- list()
lambdatilde$b2v <- list()
lambdatilde$ncobias <- list()

proportionselected <- list()
proportionselected$b2v <- list()
proportionselected$ncobias <- list()

var_ay <- vector()
EICpsipound <- matrix(0, nrow=nrow(data), ncol=length(comparisons)*V)
EICnco <- matrix(0, nrow=nrow(data), ncol=length(comparisons)*V)

bias <- list()
bias_nco <- list()

bvt <- list()
for(v in 1:length(folds)){
  message(paste("Working on fold", v, "at", Sys.time(), sep=" "))

  lambdatilde$b2v[[v]] <- vector()
  proportionselected$b2v[[v]] <- vector()
  if(is.null(NCO)==FALSE){
    lambdatilde$ncobias[[v]] <- vector()
    proportionselected$ncobias[[v]] <- vector()
  }

  #define training set
  train <- data[sort(folds[[v]]$training_set),]

  selector[[v]] <- apply_selector_func(txinrwd, train, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons, bounds)

  if(txinrwd==TRUE){
    bvt[[v]] <- bvt_txinrwd(v, selector, NCO, comparisons, train, data, fluctuation, family)
  } else {
    bvt[[v]] <- bvt_notxinrwd(v, selector, NCO, comparisons, train, data, fluctuation, family)
  }

  lambdatilde$b2v[[v]] <- comparisons[[which(bvt[[v]]$b2v==min(bvt[[v]]$b2v))]]
  proportionselected$b2v[[v]] <- which(bvt[[v]]$b2v==min(bvt[[v]]$b2v))

  if(is.null(NCO)==FALSE){
    lambdatilde$ncobias[[v]] <- comparisons[[which(bvt[[v]]$addncobias==min(bvt[[v]]$addncobias))]]
    proportionselected$ncobias[[v]] <- which(bvt[[v]]$addncobias==min(bvt[[v]]$addncobias))
  }

  for(s in 1:length(comparisons)){
    EICpsipound[,(length(comparisons)*(v-1)+s)] <- bvt[[v]]$EICpsipound[,s]
    EICnco[,(length(comparisons)*(v-1)+s)] <- bvt[[v]]$EICnco[,s]
    var_ay[(length(comparisons)*(v-1)+s)] <- bvt[[v]]$var[s]
  }


}
check <- validpreds(data, folds, V, selector, pRCT, Delta, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons)

test_that("Confirm validpreds uses rct randomization probability for s=1", {
  expect_true(all(check[[1]]$gHat1W[which(data$S==1)]==pRCT))
})

test_that("Confirm predicting missingness mechanism if missing outcomes present", {
  expect_true(all(check[[2]]$dbarAW==1)==FALSE)
})

test_that("Qbar1W=QbarAW for obs with A=1 and Qbar0W=QbarAW for obs with A=0", {
  expect_equal(check[[2]]$QbarAW[which(data$A==1)], check[[2]]$Qbar1W[which(data$A==1)])
  expect_equal(check[[2]]$QbarAW[which(data$A==0)], check[[2]]$Qbar0W[which(data$A==0)])
})

#tests for limitdistvar function
valid_initial <- check

results$ATE <- list()
results$ATE$b2v <- vector()
results$ATE$ncobias <- vector()

results$foldATEs <- list()
results$foldATEs$b2v <- vector()
results$foldATEs$ncobias <- vector()

limitdist <- limitdistvar(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons, bounds)

test_that("Solved EICs", {
  expect_equal(mean(limitdist$EICay), 0)
})

#tests for limitdist_sample
limitdistsamp <- limitdist_sample(V, bvt, NCO, EICpsipound, EICnco, var_ay, limitdist, data, comparisons)

test_that("Confirm 0 covariance between training set bias EICs and estimation set ate EICs", {
  expect_equal(c(limitdistsamp$covMat_poundplusphi[1,21],limitdistsamp$covMat_poundplusphi[1,22],limitdistsamp$covMat_poundplusphi[2,21],limitdistsamp$covMat_poundplusphi[2,22]),c(0,0,0,0))
  expect_equal(c(limitdistsamp$covMat[1,21],limitdistsamp$covMat[1,22],limitdistsamp$covMat[2,21],limitdistsamp$covMat[2,22]),c(0,0,0,0))
})

test_that("Sampling distribution approx. mean 0", {
  expect_true(abs(mean(limitdistsamp$psi_pstarnv_b2v))<0.01)
})


