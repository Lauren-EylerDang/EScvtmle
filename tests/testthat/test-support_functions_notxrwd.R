#generate data for testing
n1 <- 200
n2 <- 400

bA <- -0.7

generateData <- function(n1, n2, pRCT){
  N <- sum(n1, n2)

  #S
  S <- rep(NA, N)
  S[1:n1] <- 1
  S[(n1+1):(n1+n2)] <- 2

  #Some Ws
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)

  #A
  A <- rep(0, N)
  A[which(S==1)] <- rbinom(n1, 1, pRCT)
  A[which(S==2)] <- 0

  #Outcome regression
  Uy <- rnorm(N, 0, 1)
  Unco <- rnorm(N, 0, 1)

  Y <- 0.1*W1 - 0.5*W2 +bA*A + Uy
  nco <- 0.5*W1 - 0.1*W2 + Unco

  data <- data.frame(S, W1, W2, A, Y, nco)

  return(data)
}

data <- generateData(n1, n2, 0.5)

library(SuperLearner)

txinrwd=FALSE
study="S"
covariates=c("W1", "W2")
treatment_var="A"
treatment=1
outcome="Y"
NCO="nco"
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
bounds=NULL

data <- preprocess(txinrwd=FALSE, data=data, study="S", covariates=c("W1", "W2"), treatment_var="A", treatment=1, outcome="Y", NCO="nco", Delta=NULL, Delta_NCO=NULL, adjustnco=FALSE)
if("Delta" %in% colnames(data)){
  Delta = "Delta"
}

if("NCO_delta" %in% colnames(data)){
  Delta_NCO = "NCO_delta"
}

#Create cross-validation folds that preserve proportion of RCT (if txinwrd=TRUE) or of RCT controls (if txinrwd=FALSE) in validation sets
ids <- data$S
ids[which(data$S==1 & data$A==0)] <- 0

folds <- make_folds(data, fold_fun = folds_vfold, V=V, strata_ids = ids)
data$v <- rep(NA, nrow(data))
for(v in 1:V){
  data$v[folds[[v]]$validation_set]<-v
}

v<-1

data[which(data$S==2),]$S <- 0
train <- data[sort(folds[[v]]$training_set),]

check <- selector_func_notxrwd(train_s=train, data, Q.SL.library=c("SL.glm"), d.SL.library=NULL, g.SL.library=c("SL.glm"),
                             pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                             NCO="nco", Delta=NULL, Delta_NCO = NULL,
                             adjustnco=FALSE, target.gwt=TRUE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, bounds)


#tests for selector_func_notxrwd

test_that("Y scaled appropriately if continuous outcome with logistic fluctuation", {
  expect_true(all(check$Y>=0))
  expect_true(all(check$Y<=1))
})

test_that("Confirm QbarAW != QbarSAW", {
  expect_true(all((check$QbarAW==check$QbarSAW)==FALSE))
})

test_that("Qbar1W=QbarAW for obs with A=1 and Qbar0W=QbarAW for obs with A=0", {
  expect_equal(check$QbarAW[which(train$A==1)], check$Qbar1W[which(train$A==1)])
  expect_equal(check$QbarAW[which(train$A==0)], check$Qbar0W[which(train$A==0)])
})

test_that("Estimating missingness mechanism if missing outcomes", {
  expect_true(is.null(check$DbarSL))
  dat <- train
  dat$Y[1:10] <- NA
  dat$Delta <- c(rep(0,10),rep(1, (nrow(dat)-10)))
  out <- selector_func_notxrwd(train_s=dat, data, Q.SL.library=c("SL.glm"), d.SL.library=c("SL.glm"), g.SL.library=c("SL.glm"),
                             pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                             NCO=NULL, Delta="Delta", Delta_NCO = NULL,
                             adjustnco=FALSE, target.gwt=TRUE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, bounds)
  expect_true(is.null(out$DbarSL)==FALSE)

})

test_that("Known randomization probability used for RCT only", {
  dat1 <- data[which(data$S==1),]
  out <- selector_func_notxrwd(train_s=dat1, data, Q.SL.library=c("SL.glm"), d.SL.library=c("SL.glm"), g.SL.library=c("SL.glm"),
                             pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                             NCO=NULL, Delta=NULL, Delta_NCO = NULL,
                             adjustnco=FALSE, target.gwt=TRUE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, bounds)
  expect_equal(out$wt, rep(2, length(out$wt)))
  out <- selector_func_notxrwd(train_s=dat1, data, Q.SL.library=c("SL.glm"), d.SL.library=c("SL.glm"), g.SL.library=c("SL.glm"),
                             pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                             NCO=NULL, Delta=NULL, Delta_NCO = NULL,
                             adjustnco=FALSE, target.gwt=FALSE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, bounds)
  expect_equal(out$H.AW[which(dat1$A==1)], rep(2, length(which(dat1$A==1))))
  expect_equal(out$H.AW[which(dat1$A==0)], rep(-2, length(which(dat1$A==0))))
})

test_that("nco scaled appropriately if continuous outcome with logistic fluctuation", {
  out <- selector_func_notxrwd(train_s=train, data, Q.SL.library=c("SL.glm"), d.SL.library=NULL, g.SL.library=c("SL.glm"),
                             pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                             NCO=NULL, Delta=NULL, Delta_NCO = NULL,
                             adjustnco=FALSE, target.gwt=TRUE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, bounds)
  expect_true(all(out$nco>=0))
  expect_true(all(out$nco<=1))
})

test_that("Predictions for missingness mechanism working as expected", {
  dat1 <- train
  dat1$Delta <- c(rep(0,100), rep(1,(nrow(train)-100)))
  out <- selector_func_notxrwd(train_s=dat1, data, Q.SL.library=c("SL.glm"), d.SL.library=c("SL.glm"), g.SL.library=c("SL.glm"),
                             pRCT = 0.5, family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                             NCO=NULL, Delta="Delta", Delta_NCO = NULL,
                             adjustnco=FALSE, target.gwt=FALSE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, bounds)
  expect_equal(mean(predict(out$DbarSL, newdata = dat1[which(dat1$A==1),])), mean(dat1$Delta[which(dat1$A==1)]))
})


#tests for bvt_notxinrwd
data[which(data$S==0),]$S <- 2
train <- data[sort(folds[[v]]$training_set),]

out <- list()
out[[1]] <- apply_selector_func(txinrwd=FALSE, train=train, data, Q.SL.library=c("SL.glm"), d.SL.library=c("SL.glm"), g.SL.library=c("SL.glm"),
                                pRCT=0.5, family="gaussian", family_nco="gaussian", fluctuation="logistic", NCO="nco", Delta=NULL, Delta_NCO=NULL,
                                adjustnco=FALSE, target.gwt=TRUE, Q.discreteSL=TRUE, d.discreteSL=TRUE, g.discreteSL=TRUE, comparisons=list(c(1),c(1,2)), bounds)
check <- bvt_notxinrwd(v=1, selector=out, NCO="nco", comparisons=list(c(1),c(1,2)), train=train, data=data, fluctuation="logistic", family="gaussian")

test_that("Solves EICs", {
  expect_equal(mean(check$EIClambdav[[1]]),0)
  expect_equal(mean(check$EIClambdav[[2]]),0)
  expect_equal(mean(check$EICpsipound[,1]),0)
  expect_equal(mean(check$EICpsipound[,2]),0)
  expect_equal(mean(check$EICnco[,1]),0)
  expect_equal(mean(check$EICnco[,2]),0)
})

test_that("Bias 0 for RCT only", {
  expect_equal(check$bias[1],0)
})
