#load packages
library(tmle)
library(SuperLearner)
library(origami)
library(dplyr)
library(MASS)

#bound clever covariates
.bound <- function(x,n){
  x <- pmax((5/(n)^(1/2)/log(n)), pmin(1,x))
  return(x)
}

#Function for selecting among RCT only +/- multiple potential Real World Datasets Based on the Bias-Variance Tradeoff
selector_func <- function(train_s, data, Q.SL.library, d.SL.library, g.SL.library, pRCT = pRCT, family, family_nco, fluctuation = "logistic", NCO=NULL, Delta=NULL, Delta_NCO = NULL, adjustnco=adjustnco, target.gwt=target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL){

  #Estimate bias
  if(any(train_s$S==0)){

    Y <- train_s$Y
    if(family=="gaussian" & fluctuation == "logistic"){
      Y <- (Y - min(data$Y, na.rm = TRUE))/(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE))
    }

    if(adjustnco == FALSE){
      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "nco", "NCO_delta", "v"))==FALSE)]
    } else {
      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "NCO_delta", "v"))==FALSE)]
    }


    # set the A=0 in X0
    X0 <- X1 <- X
    X0$A <- 0
    X1$A <- 1

    if(is.null(Delta)==FALSE){
      Ynomiss <- Y[which(X$Delta==1)]
      Xnomiss <- X[which(X$Delta==1),]
      Xnomiss <- subset(Xnomiss, select=-c(Delta))
    } else {
      Ynomiss <- Y
      Xnomiss <- X
    }

    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL<- SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL$cvRisk)
        QbarSL <- QbarSL$fitLibrary[[keepAlg]]
      }

    } else {
      QbarSL<- SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family=family)
      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL$cvRisk)
        QbarSL <- QbarSL$fitLibrary[[keepAlg]]
      }
    }

    # initial estimates of the outcome, given the observed exposure & covariates
    if(Q.discreteSL==TRUE){
      QbarAW <- predict(QbarSL, newdata=X)
      # estimates of the outcome, given A=0 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)
      Qbar1W<- predict(QbarSL, newdata=X1)
    } else {
      QbarAW <- predict(QbarSL, newdata=X)$pred
      # estimates of the outcome, given A=0 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)$pred
      Qbar1W<- predict(QbarSL, newdata=X1)$pred
    }

    dHat1W <- list()
    if(is.null(Delta)==FALSE){
      DbarSL<- SuperLearner(Y=X$Delta, X=subset(X, select=-c(Delta)), SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
      if(d.discreteSL==TRUE){
        keepAlg <- which.min(DbarSL$cvRisk)
        DbarSL <- DbarSL$fitLibrary[[keepAlg]]
        dHat1W$A1 <- predict(DbarSL, newdata = X1)
        dHat1W$A0 <- predict(DbarSL, newdata = X0)
      } else {
        dHat1W$A1 <- predict(DbarSL, newdata = X1)$pred
        dHat1W$A0 <- predict(DbarSL, newdata = X0)$pred
      }
    } else {
      X$Delta <- dHat1W$A1 <- dHat1W$A0 <- rep(1, nrow(X))
      DbarSL <- NULL
    }

    # Estimate the exposure mechanism g(A|W)
    gHatSL<- SuperLearner(Y=X$A, X=subset(X, select= -c(A,Delta)), SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
    if(g.discreteSL==TRUE){
      keepAlg <- which.min(gHatSL$cvRisk)
      gHatSL <- gHatSL$fitLibrary[[keepAlg]]
      gHat1W<- predict(gHatSL, newdata = X)
    } else {
      gHat1W<- gHatSL$SL.predict
    }

    # predicted prob of not being exposed, given baseline covariates
    gHat0W<- 1- gHat1W

    #-------------------------------------------------
    # Clever covariate H(A,W) for each subject
    #-------------------------------------------------

    if(target.gwt){
      wt <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s)) + 1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s))
      H.AW1 <- as.numeric(X$A==1 & X$Delta==1)
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))
      H.AW <- (as.numeric(X$A==1 & X$Delta==1))-1*(as.numeric(X$A==0 & X$Delta==1))

      # also want to evaluate the clever covariates at A=0 for all subjects
      H.0W<- rep(-1, nrow(X))
      H.1W<- rep(1, nrow(X))
    } else{
      wt <- rep(1, nrow(X))
      H.AW1 <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s))
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s))

      H.AW <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s))-1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s))

      # also want to evaluate the clever covariates at A=0 for all subjects
      H.0W<- -1/.bound((gHat0W*dHat1W$A0),nrow(train_s))
      H.1W<- 1/.bound((gHat1W*dHat1W$A1),nrow(train_s))
    }

    #TMLE for E[E[Y|W,A=0,S=1]]
    if(adjustnco == FALSE){
      XS <- train_s[,which((colnames(train_s) %in% c("Y", "nco", "NCO_delta", "v"))==FALSE)]
    } else {
      XS <- train_s[,which((colnames(train_s) %in% c("Y", "NCO_delta", "v"))==FALSE)]
    }


    # set S=1 and A to 0 or 1
    XS1 <- XS0 <- XS
    XS1$S <- XS0$S <- 1
    XS1$A <- 1
    XS0$A <- 0

    if(is.null(Delta)==FALSE){
      XSnomiss <- XS[which(XS$Delta==1),]
      XSnomiss <- subset(XSnomiss, select=-c(Delta))
    } else {
      XSnomiss <- XS
      XS$Delta <- rep(1, nrow(XS))
    }
    # call Super Learner for estimation of QbarSAW
    if(fluctuation == "logistic"){
      QbarSL_S<- SuperLearner(Y=Ynomiss, X=XSnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))

      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL_S$cvRisk)
        QbarSL_S <- QbarSL_S$fitLibrary[[keepAlg]]
      }
    } else {
      QbarSL_S<- SuperLearner(Y=Ynomiss, X=XSnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE))

      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL_S$cvRisk)
        QbarSL_S <- QbarSL_S$fitLibrary[[keepAlg]]
      }
    }

    # initial estimates of the outcome, given the observed exposure & covariates
    if(Q.discreteSL==TRUE){
      QbarSAW <- predict(QbarSL_S, newdata=XS)
      # estimates of the outcome, given S=1 and A=1 or A=0 and covariates
      QbarS1W<- predict(QbarSL_S, newdata=XS1)
      QbarS0W<- predict(QbarSL_S, newdata=XS0)
    } else {
      QbarSAW <- predict(QbarSL_S, newdata=XS)$pred
      # estimates of the outcome, given S=1 and A=1 or A=0 and covariates
      QbarS1W<- predict(QbarSL_S, newdata=XS1)$pred
      QbarS0W<- predict(QbarSL_S, newdata=XS0)$pred
    }

    # Estimate the trial participation mechanism g(S|A=0,Delta=1,W)
    #------------------------------------------
    XS_D1 <- XS[which(XS$Delta==1),]
    # call Super Learner for the exposure mechanism
    gSHatSL<- SuperLearner(Y=XS_D1$S, X=subset(XS_D1, select= -c(S,Delta)), SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
    if(g.discreteSL==TRUE){
      keepAlg <- which.min(gSHatSL$cvRisk)
      gSHatSL <- gSHatSL$fitLibrary[[keepAlg]]
      gSHat1W <- predict(gSHatSL, newdata = XS1)
      gSHat0W <- predict(gSHatSL, newdata = XS0)
    } else {
      # generate predicted prob being exposed, given baseline covariates
      gSHat1W <- predict(gSHatSL, newdata = XS1)$pred
      gSHat0W <- predict(gSHatSL, newdata = XS0)$pred
    }


    #-------------------------------------------------
    # Clever covariate H(S,A,W) for each subject
    #-------------------------------------------------
    if(target.gwt){
      wt_s <- as.numeric(XS$A==1 & XS$S==1 & XS$Delta == 1)/.bound((gSHat1W*gHat1W*dHat1W$A1),nrow(train_s)) + as.numeric(XS$A==0 & XS$S==1 & XS$Delta == 1)/.bound((gSHat0W*gHat0W*dHat1W$A0),nrow(train_s))
      H.SAW1 <- as.numeric(XS$A==1 & XS$S==1 & XS$Delta == 1)
      H.SAW0 <- -1*as.numeric(XS$A==0 & XS$S==1 & XS$Delta == 1)

      # also want to evaluate the clever covariates at S=1 and A=1 or 0 for all subjects
      H.S1W<- rep(1, nrow(XS))
      H.S0W<- rep(-1, nrow(XS))
    } else{
      wt_s <- rep(1, nrow(XS))
      H.SAW1 <- as.numeric(XS$A==1 & XS$S==1 & XS$Delta == 1)/.bound((gSHat1W*gHat1W*dHat1W$A1),nrow(train_s))
      H.SAW0 <- -1*as.numeric(XS$A==0 & XS$S==1 & XS$Delta == 1)/.bound((gSHat0W*gHat0W*dHat1W$A0),nrow(train_s))

      # also want to evaluate the clever covariates at S=1 and A=1 or 0 for all subjects
      H.S1W<- 1/.bound((gSHat1W*gHat1W*dHat1W$A1),nrow(train_s))
      H.S0W<- -1/.bound((gSHat0W*gHat0W*dHat1W$A0),nrow(train_s))
    }

    if(is.null(NCO)==FALSE){
      # call Super Learner for estimation of NCObarAW
      if(family_nco=="gaussian" & fluctuation == "logistic"){
        train_s_nco <- (train_s$nco - min(data$nco, na.rm=TRUE))/(max(data$nco, na.rm=TRUE) - min(data$nco, na.rm=TRUE))
      } else {
        train_s_nco <- train_s$nco
      }

      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "nco", "Delta", "v"))==FALSE)]
      X0 <- X1 <- X
      X0$A <- 0
      X1$A <- 1

      if(is.null(Delta_NCO)==FALSE){
        NCOnomiss <- train_s_nco[which(X$NCO_delta==1)]
        Xnomiss <- X[which(X$NCO_delta==1),]
        Xnomiss <- subset(Xnomiss, select=-c(NCO_delta))
      } else {
        NCOnomiss <- train_s_nco
        Xnomiss <- X
        X$NCO_delta <- rep(1, nrow(X))
      }

      if(fluctuation == "logistic"){
        NCObarSL<- SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
        if(Q.discreteSL==TRUE){
          keepAlg <- which.min(NCObarSL$cvRisk)
          NCObarSL <- NCObarSL$fitLibrary[[keepAlg]]
        }
      } else {
        NCObarSL<- SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE))
        if(Q.discreteSL==TRUE){
          keepAlg <- which.min(NCObarSL$cvRisk)
          NCObarSL <- NCObarSL$fitLibrary[[keepAlg]]
        }
      }

      if(Q.discreteSL==TRUE){
        NCObarAW <- predict(NCObarSL, X)
        NCObar1W <- predict(NCObarSL, X1)
        NCObar0W <- predict(NCObarSL, X0)
      } else {
        NCObarAW <- predict(NCObarSL, X)$pred
        NCObar1W <- predict(NCObarSL, X1)$pred
        NCObar0W <- predict(NCObarSL, X0)$pred
      }

      # Estimate the exposure mechanism g(A|W)
      gHatSLnco<- SuperLearner(Y=X$A, X=subset(X, select= -c(A,NCO_delta)), SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
      if(g.discreteSL==TRUE){
        keepAlg <- which.min(gHatSLnco$cvRisk)
        gHatSLnco <- gHatSLnco$fitLibrary[[keepAlg]]
        gHat1W<- predict(gHatSLnco, X)
      } else {
        gHat1W<- gHatSLnco$SL.predict
      }

      # predicted prob of not being exposed, given baseline covariates
      gHat0W<- 1- gHat1W

      dHat1Wnco <- list()
      if(is.null(Delta_NCO)==FALSE){
        DbarSLnco<- SuperLearner(Y=X$NCO_delta, X=subset(X, select=-c(NCO_delta)), SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
        if(d.discreteSL==TRUE){
          keepAlg <- which.min(DbarSLnco$cvRisk)
          DbarSLnco <- DbarSLnco$fitLibrary[[keepAlg]]
          dHat1Wnco$A1 <- predict(DbarSLnco, newdata = X1)
          dHat1Wnco$A0 <- predict(DbarSLnco, newdata = X0)
        } else {
          dHat1Wnco$A1 <- predict(DbarSLnco, newdata = X1)$pred
          dHat1Wnco$A0 <- predict(DbarSLnco, newdata = X0)$pred
        }
      } else {
        dHat1Wnco$A1 <- dHat1Wnco$A0 <- rep(1, nrow(X))
        DbarSLnco <- NULL
      }

      #-------------------------------------------------
      # Clever covariate H(A,W) for each subject
      #-------------------------------------------------
      if(target.gwt){
        wt_nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s)) + 1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))
        H.AW1nco <- as.numeric(X$A==1 & X$NCO_delta==1)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))
        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))-1*(as.numeric(X$A==0 & X$NCO_delta==1))

        # also want to evaluate the clever covariates at A=0 for all subjects
        H.0Wnco<- rep(-1, nrow(X))
        H.1Wnco<- rep(1, nrow(X))
      } else{
        wt_nco <- rep(1, nrow(X))
        H.AW1nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s))
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))

        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s))-1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))

        # also want to evaluate the clever covariates at A=0 for all subjects
        H.0Wnco<- -1/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))
        H.1Wnco<- 1/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s))
      }

    } else {
      NCObarAW <- NULL
      NCObar1W <- NULL
      NCObar0W <- NULL
      H.AWnco <- NULL
      H.0Wnco <- NULL
      H.1Wnco <- NULL
      train_s_nco <- NULL
      wt_nco <- NULL
    }

  } else {
    Y <- train_s$Y
    if(family=="gaussian" & fluctuation == "logistic"){
      Y <- (Y - min(data$Y, na.rm = TRUE))/(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE))
    }

    #run tmle for E[E[Y|W,A=0]]
    if(adjustnco == FALSE){
      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "nco", "NCO_delta", "v"))==FALSE)]
    } else {
      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "NCO_delta", "v"))==FALSE)]
    }


    # set the A=0 in X0
    X0 <- X1 <- X
    X0$A <- 0
    X1$A <- 1

    if(is.null(Delta)==FALSE){
      Ynomiss <- Y[which(X$Delta==1)]
      Xnomiss <- X[which(X$Delta==1),]
      Xnomiss <- subset(Xnomiss, select=-c(Delta))
    } else {
      Ynomiss <- Y
      Xnomiss <- X
    }

    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL<- SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL$cvRisk)
        QbarSL <- QbarSL$fitLibrary[[keepAlg]]
      }
    } else {
      QbarSL<- SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE))
      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL$cvRisk)
        QbarSL <- QbarSL$fitLibrary[[keepAlg]]
      }
    }

    if(Q.discreteSL==TRUE){
      # initial estimates of the outcome, given the observed exposure & covariates
      QbarAW <- predict(QbarSL, newdata=X)
      # estimates of the outcome, given A=0 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)
      Qbar1W<- predict(QbarSL, newdata=X1)
    } else {
      # initial estimates of the outcome, given the observed exposure & covariates
      QbarAW <- predict(QbarSL, newdata=X)$pred
      # estimates of the outcome, given A=0 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)$pred
      Qbar1W<- predict(QbarSL, newdata=X1)$pred
    }

    dHat1W <- list()
    if(is.null(Delta)==FALSE){
      DbarSL<- SuperLearner(Y=X$Delta, X=subset(X, select=-c(Delta)), SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
      if(d.discreteSL==TRUE){
        keepAlg <- which.min(DbarSL$cvRisk)
        DbarSL <- DbarSL$fitLibrary[[keepAlg]]
        dHat1W$A1 <- predict(DbarSL, newdata = X1)
        dHat1W$A0 <- predict(DbarSL, newdata = X0)
      } else {
        dHat1W$A1 <- predict(DbarSL, newdata = X1)$pred
        dHat1W$A0 <- predict(DbarSL, newdata = X0)$pred
      }
    } else {
      X$Delta <- dHat1W$A1 <- dHat1W$A0 <- rep(1, nrow(X))
      DbarSL <- NULL
    }

    # Estimate the exposure mechanism g(A|W)
    gHatSL<- NULL
    gHat1W<- rep(pRCT, nrow(X))
    # predicted prob of not being exposed, given baseline covariates
    gHat0W<- 1- gHat1W

    #-------------------------------------------------
    # Clever covariate H(A,W) for each subject
    #-------------------------------------------------
    if(target.gwt){
      wt <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s)) + 1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s))
      H.AW1 <- as.numeric(X$A==1 & X$Delta==1)
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))
      H.AW <- (as.numeric(X$A==1 & X$Delta==1))-1*(as.numeric(X$A==0 & X$Delta==1))

      # also want to evaluate the clever covariates at A=0 and A=1 for all subjects
      H.0W<- rep(-1, nrow(X))
      H.1W<- rep(1, nrow(X))
    } else{
      wt <- rep(1, nrow(X))
      H.AW1 <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s))
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s))

      H.AW <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s))-1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s))

      # also want to evaluate the clever covariates at A=0 for all subjects
      H.0W<- -1/.bound((gHat0W*dHat1W$A0),nrow(train_s))
      H.1W<- 1/.bound((gHat1W*dHat1W$A1),nrow(train_s))
    }

    QbarSAW <- NULL
    QbarS1W <- NULL
    QbarS0W <- NULL
    H.SAW1 <- NULL
    H.SAW0 <- NULL
    H.S1W <- NULL
    H.S0W <- NULL
    wt_s <- rep(1, nrow(X))

    if(is.null(NCO)==FALSE){
      # call Super Learner for estimation of NCObarAW
      if(family_nco=="gaussian" & fluctuation == "logistic"){
        train_s_nco <- (train_s$nco - min(data$nco, na.rm=TRUE))/(max(data$nco, na.rm=TRUE) - min(data$nco, na.rm=TRUE))
      } else {
        train_s_nco <- train_s$nco
      }

      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "nco", "Delta", "v"))==FALSE)]
      X0 <- X1 <- X
      X0$A <- 0
      X1$A <- 1

      if(is.null(Delta_NCO)==FALSE){
        NCOnomiss <- train_s_nco[which(X$NCO_delta==1)]
        Xnomiss <- X[which(X$NCO_delta==1),]
        Xnomiss <- subset(Xnomiss, select=-c(NCO_delta))
      } else {
        NCOnomiss <- train_s_nco
        Xnomiss <- X
      }

      if(fluctuation == "logistic"){
        NCObarSL<- SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
        if(Q.discreteSL==TRUE){
          keepAlg <- which.min(NCObarSL$cvRisk)
          NCObarSL <- NCObarSL$fitLibrary[[keepAlg]]
        }
      } else {
        NCObarSL<- SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE))
        if(Q.discreteSL==TRUE){
          keepAlg <- which.min(NCObarSL$cvRisk)
          NCObarSL <- NCObarSL$fitLibrary[[keepAlg]]
        }
      }

      if(Q.discreteSL==TRUE){
        NCObarAW <- predict(NCObarSL, X)
        NCObar1W <- predict(NCObarSL, X1)
        NCObar0W <- predict(NCObarSL, X0)
      } else {
        NCObarAW <- predict(NCObarSL, X)$pred
        NCObar1W <- predict(NCObarSL, X1)$pred
        NCObar0W <- predict(NCObarSL, X0)$pred
      }




      gHat1W<- rep(pRCT, nrow(X))
      # predicted prob of not being exposed, given baseline covariates
      gHat0W<- 1- gHat1W

      dHat1Wnco <- list()
      if(is.null(Delta_NCO)==FALSE){
        DbarSLnco<- SuperLearner(Y=X$NCO_delta, X=subset(X, select=-c(NCO_delta)), SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
        if(d.discreteSL==TRUE){
          keepAlg <- which.min(DbarSLnco$cvRisk)
          DbarSLnco <- DbarSLnco$fitLibrary[[keepAlg]]
          dHat1Wnco$A1 <- predict(DbarSLnco, newdata = X1)
          dHat1Wnco$A0 <- predict(DbarSLnco, newdata = X0)
        } else {
          dHat1Wnco$A1 <- predict(DbarSLnco, newdata = X1)$pred
          dHat1Wnco$A0 <- predict(DbarSLnco, newdata = X0)$pred
        }
      } else {
        X$NCO_delta <- dHat1Wnco$A1 <- dHat1Wnco$A0 <- rep(1, nrow(X))
        DbarSLnco <- NULL
      }

      #-------------------------------------------------
      # Clever covariate H(A,W) for each subject
      #-------------------------------------------------
      if(target.gwt){
        wt_nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s)) + 1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))
        H.AW1nco <- as.numeric(X$A==1 & X$NCO_delta==1)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))
        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))-1*(as.numeric(X$A==0 & X$NCO_delta==1))

        # also want to evaluate the clever covariates at A=0 for all subjects
        H.0Wnco<- rep(-1, nrow(X))
        H.1Wnco<- rep(1, nrow(X))
      } else{
        wt_nco <- rep(1, nrow(X))
        H.AW1nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s))
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))

        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s))-1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))

        # also want to evaluate the clever covariates at A=0 for all subjects
        H.0Wnco<- -1/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s))
        H.1Wnco<- 1/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s))
      }

    } else {
      NCObarAW <- NULL
      NCObar1W <- NULL
      NCObar0W <- NULL
      H.AWnco <- NULL
      H.0Wnco <- NULL
      H.1Wnco <- NULL
      train_s_nco <- NULL
      wt_nco <- NULL
    }
  }

  out <- list("Y" = Y, "QbarSL" = QbarSL, "QbarAW" = QbarAW, "Qbar1W" = Qbar1W, "Qbar0W" = Qbar0W, "QbarSAW" = QbarSAW, "QbarS1W" = QbarS1W, "QbarS0W" = QbarS0W, "gHatSL" = gHatSL, "train_s_nco" = train_s_nco, "NCObarAW" = NCObarAW, "NCObar1W" = NCObar1W, "NCObar0W" = NCObar0W, "H.AW1" = H.AW1, "H.AW0" = H.AW0,"H.AW" = H.AW, "H.0W" = H.0W, "H.1W" = H.1W, "H.SAW1" = H.SAW1,"H.SAW0" = H.SAW0, "H.S1W" = H.S1W, "H.S0W" = H.S0W, "DbarSL" = DbarSL, "H.AWnco" = H.AWnco, "H.0Wnco" = H.0Wnco, "H.1Wnco"= H.1Wnco, "wt"=wt, "wt_s"=wt_s, "wt_nco"=wt_nco)

  return(out)

}

#apply selector_func to different datasets
apply_selector_func <- function(train, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO=NULL, Delta=NULL, Delta_NCO=NULL, adjustnco=adjustnco, target.gwt=target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL, comparisons){
  out <- list()
  for(s in 1:(length(comparisons))){

    #train regressions on A=0 only
    train_s <- train[which(train$S %in% comparisons[[s]]),]
    train_s$S[which(train_s$S!=1)]<-0

    out[[s]] <- selector_func(train_s = train_s, data=data, Q.SL.library, d.SL.library, g.SL.library, pRCT = pRCT, family = family, family_nco = family_nco, fluctuation = fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL)
  }
  return(out)
}



validpreds <- function(data, folds, V, selector, pRCT, Delta=NULL, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons){
  out <- list()
  for(s in 1:length(comparisons)){
    out[[s]] <- matrix(0, nrow=nrow(data), ncol=8)
    out[[s]] <- data.frame(out[[s]])
    colnames(out[[s]]) <- c("v", "QbarAW", "Qbar1W", "Qbar0W", "dbarAW", "dbar1W", "dbar0W", "gHat1W")

    out[[s]]$v <- data$v
    out[[s]]$dbarAW <- out[[s]]$dbar1W <- out[[s]]$dbar0W <- rep(1, nrow(out[[s]]))

    D1 <- D0 <- data
    D1$A <- 1
    D0$A <- 0

    for(v in 1:V){
      if(Q.discreteSL==TRUE){
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$QbarAW <- predict(selector[[v]][[s]]$QbarSL, newdata = subset(data[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y)))
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar1W <- predict(selector[[v]][[s]]$QbarSL, newdata = subset(D1[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y)))
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar0W <- predict(selector[[v]][[s]]$QbarSL, newdata = subset(D0[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y)))
      } else {
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$QbarAW <- predict(selector[[v]][[s]]$QbarSL, newdata = subset(data[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y)))$pred
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar1W <- predict(selector[[v]][[s]]$QbarSL, newdata = subset(D1[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y)))$pred
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar0W <- predict(selector[[v]][[s]]$QbarSL, newdata = subset(D0[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y)))$pred
      }

      if(is.null(Delta)==FALSE){
        if(d.discreteSL==TRUE){
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbarAW <- predict(selector[[v]][[s]]$DbarSL, newdata = subset(data[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y,Delta)))
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar1W <- predict(selector[[v]][[s]]$DbarSL, newdata = subset(D1[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y,Delta)))
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar0W <- predict(selector[[v]][[s]]$DbarSL, newdata = subset(D0[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y,Delta)))
        } else {
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbarAW <- predict(selector[[v]][[s]]$DbarSL, newdata = subset(data[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y,Delta)))$pred
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar1W <- predict(selector[[v]][[s]]$DbarSL, newdata = subset(D1[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y,Delta)))$pred
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar0W <- predict(selector[[v]][[s]]$DbarSL, newdata = subset(D0[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(S,Y,Delta)))$pred
        }
      }

      if(s==1){
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$gHat1W <- rep(pRCT, length(which(data$v==v & (data$S %in% comparisons[[s]]))))
      } else {
        if(length(comparisons)>1){
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$gHat1W <- predict(selector[[v]][[s]]$gHatSL, newdata = subset(data[which(data$v==v & (data$S %in% comparisons[[s]])),], select=-c(A,S,Y)))
        }
      }
    }
  }
  return(out)
}

#function for limit dist var ests
limitdistvar<- function(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons){
  out <- list()

  out$EICay <- matrix(0, nrow=nrow(data), ncol=length(comparisons)*V)
  out$psi <- list()
  for(v in 1:V){
    out$psi[[v]] <- vector()
  }

  out$psi_s <- vector()

  for(s in 1:length(comparisons)){

    #select experiment subset

    validmat <- valid_initial[[s]]
    validmat$Yscale <- data$Y

    if(family=="gaussian" & fluctuation == "logistic"){
      validmat$Yscale <- (validmat$Yscale - min(data$Y))/(max(data$Y) - min(data$Y))
    }

    if(all(is.na(str_match(colnames(data), "Delta"))==TRUE)){
      data$Delta <- rep(1, nrow(data))
    }

    if(target.gwt){
      wt <- as.numeric(data$A==1 & data$Delta==1)/.bound((validmat$gHat1W*validmat$dbar1W),length(which(data$S %in% comparisons[[s]]))) + as.numeric(data$A==0 & data$Delta==1)/.bound(((1-validmat$gHat1W)*validmat$dbar0W),length(which(data$S %in% comparisons[[s]])))
      H.AW <- as.numeric(data$A==1 & data$Delta==1) - as.numeric(data$A==0 & data$Delta==1)
      H.1W <- rep(1, nrow(data))
      H.0W <- rep(-1, nrow(data))

    } else{
      wt <- rep(1, nrow(data))
      H.AW <- as.numeric(data$A==1 & data$Delta==1)/.bound((validmat$gHat1W*validmat$dbar1W),length(which(data$S %in% comparisons[[s]]))) - as.numeric(data$A==0 & data$Delta==1)/.bound(((1-validmat$gHat1W)*validmat$dbar0W),length(which(data$S %in% comparisons[[s]])))
      H.1W <- 1/.bound((validmat$gHat1W*validmat$dbar1W),length(which(data$S %in% comparisons[[s]])))
      H.0W <- -1/.bound(((1-validmat$gHat1W)*validmat$dbar0W),length(which(data$S %in% comparisons[[s]])))

    }


    if(fluctuation == "logistic"){
      logitUpdate<- glm(validmat$Yscale[which(data$Delta==1 & (data$S %in% comparisons[[s]]))] ~ -1 + offset(qlogis(validmat$QbarAW[which(data$Delta==1 & (data$S %in% comparisons[[s]]))])) +  H.AW[which(data$Delta==1 & (data$S %in% comparisons[[s]]))], family='binomial', weights = wt[which(data$Delta==1 & (data$S %in% comparisons[[s]]))])

      epsilon <- logitUpdate$coef

      validmat$QbarAW.star <- validmat$Qbar1W.star <- validmat$Qbar0W.star <- rep(0, nrow(validmat))

      validmat[which(data$S %in% comparisons[[s]]),]$QbarAW.star <- plogis(qlogis(validmat[which(data$S %in% comparisons[[s]]),]$QbarAW) + epsilon*H.AW[which(data$S %in% comparisons[[s]])])
      validmat[which(data$S %in% comparisons[[s]]),]$Qbar1W.star<- plogis(qlogis(validmat[which(data$S %in% comparisons[[s]]),]$Qbar1W) + epsilon*H.1W[which(data$S %in% comparisons[[s]])])
      validmat[which(data$S %in% comparisons[[s]]),]$Qbar0W.star<- plogis(qlogis(validmat[which(data$S %in% comparisons[[s]]),]$Qbar0W) + epsilon*H.0W[which(data$S %in% comparisons[[s]])])
      if(family == "gaussian"){
        validmat[which(data$S %in% comparisons[[s]]),]$QbarAW.star <- validmat[which(data$S %in% comparisons[[s]]),]$QbarAW.star*(max(data$Y) - min(data$Y)) + min(data$Y)
        validmat[which(data$S %in% comparisons[[s]]),]$Qbar1W.star <- validmat[which(data$S %in% comparisons[[s]]),]$Qbar1W.star*(max(data$Y) - min(data$Y)) + min(data$Y)
        validmat[which(data$S %in% comparisons[[s]]),]$Qbar0W.star <- validmat[which(data$S %in% comparisons[[s]]),]$Qbar0W.star*(max(data$Y) - min(data$Y)) + min(data$Y)
      }

    } else {
      logitUpdate<- glm(validmat$Yscale[which(data$Delta==1 & (data$S %in% comparisons[[s]]))] ~ -1 + offset(validmat$QbarAW[which(data$Delta==1 & (data$S %in% comparisons[[s]]))]) +  H.AW[which(data$Delta==1 & (data$S %in% comparisons[[s]]))], family='gaussian', weights = wt[which(data$Delta==1 & (data$S %in% comparisons[[s]]))])

      epsilon <- logitUpdate$coef

      validmat$QbarAW.star <- validmat$Qbar1W.star <- validmat$Qbar0W.star <- rep(0, nrow(validmat))

      validmat[which(data$S %in% comparisons[[s]]),]$QbarAW.star<- validmat[which(data$S %in% comparisons[[s]]),]$QbarAW + epsilon*H.AW[which(data$S %in% comparisons[[s]])]
      validmat[which(data$S %in% comparisons[[s]]),]$Qbar1W.star<- validmat[which(data$S %in% comparisons[[s]]),]$Qbar1W + epsilon*H.1W[which(data$S %in% comparisons[[s]])]
      validmat[which(data$S %in% comparisons[[s]]),]$Qbar0W.star<- validmat[which(data$S %in% comparisons[[s]]),]$Qbar0W + epsilon*H.0W[which(data$S %in% comparisons[[s]])]
    }

    for(v in 1:V){
      out$psi[[v]][s] <- mean((validmat$Qbar1W.star - validmat$Qbar0W.star)[which(validmat$v==v & (data$S %in% comparisons[[s]]))])
      out$EICay[which(data$v==v & (data$S %in% comparisons[[s]])),(length(comparisons)*(v-1)+s)] <- ((wt*H.AW)*(data$Y - validmat$QbarAW.star) + validmat$Qbar1W.star - validmat$Qbar0W.star - out$psi[[v]][s])[which(data$v==v & (data$S %in% comparisons[[s]]))]/((length(which(data$v==v & (data$S %in% comparisons[[s]]))))/nrow(data))
    }

    if(s==1){
      pooledVar <- list()
      for(v in 1:V){
        pooledVar[[v]] <- var(((wt*H.AW)*(data$Y - validmat$QbarAW.star) + validmat$Qbar1W.star - validmat$Qbar0W.star - out$psi[[v]][s])[which(data$v==v & (data$S %in% comparisons[[s]]))])/((length(which(data$S %in% comparisons[[s]]))))
      }
      out$Var <- mean(unlist(pooledVar))
    }
    out$psi_s[s] <- mean((validmat$Qbar1W.star - validmat$Qbar0W.star)[which(data$S %in% comparisons[[s]])])

  }

  return(out)
}

#add option to compare multiple different subsets of lambdas for simulation
IORD_cvtmle_sims <- function(data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, pRCT, V, Q.SL.library, d.SL.library, g.SL.library, Q.discreteSL, d.discreteSL, g.discreteSL, family, family_nco, fluctuation = "logistic", comparisons = list(c(0,1)), adjustnco = TRUE, target.gwt = TRUE){

  #remove observations missing treatment
  data <- rename(data, A = treatment_var)
  data <- data[which(is.na(data$A)==FALSE),]

  #make A coded as 1=treatment of interest, 0=control
  data$A <- ifelse(data$A == treatment, 1, 0)

  data <- rename(data, S = study)
  data <- rename(data, Y = outcome)

  #trim data to avoid positivity violation
  for(w in 1:length(covariates)){
    whichW <- which(data$S!=1 & (data[,covariates[w]] < (min(data[which(data$S==1),covariates[w]]))) | (data[,covariates[w]] > (max(data[which(data$S==1),covariates[w]]))))
    if(length(whichW)>0){
      data <- data[-whichW,]
    }
  }

  if(is.null(NCO) == FALSE){
    data <- rename(data, nco = NCO)

    if(adjustnco == TRUE){
      whichW <- which(data$S!=1 & (data[,"nco"] < (min(data[which(data$S==1),"nco"]))) | (data[,"nco"] > (max(data[which(data$S==1),"nco"]))))
      if(length(whichW)>0){
        data <- data[-whichW,]
      }
    }

    if(is.null(Delta_NCO)==FALSE){
      data <- rename(data, NCO_delta = Delta_NCO)
    }

    if(any(is.na(data$nco)==TRUE)){
      if(is.null(Delta_NCO)==TRUE){
        data$NCO_delta <- rep(1, nrow(data))
        data$NCO_delta[which(is.na(data$nco)==TRUE)] <- 0
        Delta_NCO <- "NCO_delta"
      }
      data$nco[which(is.na(data$nco)==TRUE)] <- mean(data$nco, na.rm=TRUE)
    }
  }


  if(is.null(Delta)==FALSE){
    data <- rename(data, Delta = Delta)
  }

  if(any(is.na(data$Y)==TRUE)){
    if(is.null(Delta)==TRUE){
      data$Delta <- rep(1, nrow(data))
      data$Delta[which(is.na(data$Y)==TRUE)] <- 0
      Delta <- "Delta"
    }
    data$Y[which(is.na(data$Y)==TRUE)] <- mean(data$Y, na.rm=TRUE)
  }



  data <- data[,which(colnames(data) %in% c("S", covariates, "A", "Y", "nco", "NCO_delta", "Delta"))]

  ids <- data$S
  ids[which(data$S==1)] <- 0
  folds <- make_folds(data, fold_fun = folds_vfold, V=V, strata_ids = ids)
  data$v <- rep(NA, nrow(data))
  for(v in 1:V){
    data$v[folds[[v]]$validation_set]<-v
  }

  results <- list()
  selector <- list()
  results$pooled <- list()
  valid_initial <- list()

  lambdatilde <- list()
  lambdatilde$b2v <- list()
  lambdatilde$ncobias <- list()
  lambdatilde$ncoonly <- list()

  proportionselected <- list()
  proportionselected$b2v <- list()
  proportionselected$ncobias <- list()
  proportionselected$ncoonly <- list()

  EICay <- vector()
  EICpsipound <- matrix(0, nrow=nrow(data), ncol=length(comparisons)*V)
  EICnco <- matrix(0, nrow=nrow(data), ncol=length(comparisons)*V)


  bias <- list()
  bias_nco <- list()
  for(v in 1:length(folds)){
    print(v)
    print(Sys.time())

    lambdatilde$b2v[[v]] <- vector()
    proportionselected$b2v[[v]] <- vector()
    if(is.null(NCO)==FALSE){
      lambdatilde$ncobias[[v]] <- vector()
      lambdatilde$ncoonly[[v]] <- vector()

      proportionselected$ncobias[[v]] <- vector()
      proportionselected$ncoonly[[v]] <- vector()
    }


    #define training set
    train <- data[sort(folds[[v]]$training_set),]

    selector[[v]] <- apply_selector_func(train, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons)

    b2v <- vector()
    bias[[v]] <- vector()
    var <- vector()

    if(is.null(NCO)==FALSE){
      addncobias <- vector()
      ncoonly <- vector()
      bias_nco[[v]] <- vector()
    }


    for(s in 1:length(comparisons)){
      train_s <- train[which(train$S %in% comparisons[[s]]),]
      if(is.null(train_s$Delta)){
        train_s$Delta <- rep(1, nrow(train_s))
      }

      #tmle
      if(fluctuation == "logistic"){
        logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$QbarAW[which(train_s$Delta==1)])) +  selector[[v]][[s]]$H.AW1[which(train_s$Delta==1)] + selector[[v]][[s]]$H.AW0[which(train_s$Delta==1)], family='binomial', weights = selector[[v]][[s]]$wt[which(train_s$Delta==1)])

        epsilon <- logitUpdate$coef

        QbarAW.star<- plogis(qlogis(selector[[v]][[s]]$QbarAW) + epsilon[1]*selector[[v]][[s]]$H.AW1 + epsilon[2]*selector[[v]][[s]]$H.AW0)
        Qbar1W.star<- plogis(qlogis(selector[[v]][[s]]$Qbar1W) + epsilon[1]*selector[[v]][[s]]$H.1W)
        Qbar0W.star<- plogis(qlogis(selector[[v]][[s]]$Qbar0W) + epsilon[2]*selector[[v]][[s]]$H.0W)

      } else {
        logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(selector[[v]][[s]]$QbarAW[which(train_s$Delta==1)]) +  selector[[v]][[s]]$H.AW1[which(train_s$Delta==1)] + selector[[v]][[s]]$H.AW0[which(train_s$Delta==1)], family='gaussian', weights = selector[[v]][[s]]$wt[which(train_s$Delta==1)])

        epsilon <- logitUpdate$coef

        QbarAW.star<- selector[[v]][[s]]$QbarAW + epsilon[1]*selector[[v]][[s]]$H.AW1 + epsilon[2]*selector[[v]][[s]]$H.AW0
        Qbar1W.star<- selector[[v]][[s]]$Qbar1W + epsilon[1]*selector[[v]][[s]]$H.1W
        Qbar0W.star<- selector[[v]][[s]]$Qbar0W + epsilon[2]*selector[[v]][[s]]$H.0W
      }

      if(family=="gaussian" & fluctuation == "logistic"){
        QbarAW.star <- QbarAW.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        Qbar1W.star <- Qbar1W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        Qbar0W.star <- Qbar0W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        trainsY <- selector[[v]][[s]]$Y*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
      } else {
        trainsY <- selector[[v]][[s]]$Y
      }

      EIClambdav <- (selector[[v]][[s]]$wt*(selector[[v]][[s]]$H.AW1 + selector[[v]][[s]]$H.AW0)*(trainsY - QbarAW.star) + Qbar1W.star - Qbar0W.star - mean(Qbar1W.star - Qbar0W.star))

      if(s==1){
        QbarS1W.star <- Qbar1W.star
        QbarS0W.star <- Qbar0W.star
      } else {
        #------------------------------------------
        # Update the initial estimator
        #------------------------------------------
        if(fluctuation == "logistic"){
          logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$QbarSAW[which(train_s$Delta==1)])) +  selector[[v]][[s]]$H.SAW1[which(train_s$Delta==1)] + selector[[v]][[s]]$H.SAW0[which(train_s$Delta==1)], family='binomial', weights = selector[[v]][[s]]$wt_s[which(train_s$Delta==1)])

          epsilon <- logitUpdate$coef

          QbarSAW.star<- plogis(qlogis(selector[[v]][[s]]$QbarSAW) + epsilon[1]*selector[[v]][[s]]$H.SAW1 + epsilon[2]*selector[[v]][[s]]$H.SAW0)
          QbarS1W.star<- plogis(qlogis(selector[[v]][[s]]$QbarS1W) + epsilon[1]*selector[[v]][[s]]$H.S1W)
          QbarS0W.star<- plogis(qlogis(selector[[v]][[s]]$QbarS0W) + epsilon[2]*selector[[v]][[s]]$H.S0W)
        } else {
          logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(selector[[v]][[s]]$QbarSAW[which(train_s$Delta==1)]) +  selector[[v]][[s]]$H.SAW1[which(train_s$Delta==1)] + selector[[v]][[s]]$H.SAW0[which(train_s$Delta==1)], family='gaussian', weights = selector[[v]][[s]]$wt_s[which(train_s$Delta==1)])

          epsilon <- logitUpdate$coef

          QbarSAW.star<- selector[[v]][[s]]$QbarSAW + epsilon[1]*selector[[v]][[s]]$H.SAW1 + epsilon[2]*selector[[v]][[s]]$H.SAW0
          QbarS1W.star<- selector[[v]][[s]]$QbarS1W + epsilon[1]*selector[[v]][[s]]$H.S1W
          QbarS0W.star<- selector[[v]][[s]]$QbarS0W + epsilon[2]*selector[[v]][[s]]$H.S0W
        }
        #------------------------------------------
        # Estimate Psi(P_0)
        #------------------------------------------

        if(family=="gaussian" & fluctuation == "logistic"){
          QbarSAW.star <- QbarSAW.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
          QbarS1W.star <- QbarS1W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
          QbarS0W.star <- QbarS0W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        }
      }

      if(s>1){
        EICpsipound[which(data$v!=v & data$S %in% comparisons[[s]]),(length(comparisons)*(v-1)+s)] <- (selector[[v]][[s]]$wt*(selector[[v]][[s]]$H.AW1 + selector[[v]][[s]]$H.AW0)*(trainsY - QbarAW.star) + Qbar1W.star - Qbar0W.star - mean(Qbar1W.star - Qbar0W.star) - (selector[[v]][[s]]$wt_s*(selector[[v]][[s]]$H.SAW1 + selector[[v]][[s]]$H.SAW0)*(trainsY - QbarSAW.star) + QbarS1W.star - QbarS0W.star - mean(QbarS1W.star - QbarS0W.star)))/(length(which(train$S %in% comparisons[[s]]))/nrow(data))
      }

      bias[[v]][s] <- mean(Qbar1W.star - Qbar0W.star) - mean(QbarS1W.star - QbarS0W.star)
      var[s] <- EICay[(length(comparisons)*(v-1)+s)] <- var(EIClambdav)/length(trainsY)

      #nco
      if(is.null(NCO)==FALSE){

        if(is.null(train_s$NCO_delta)){
          train_s$NCO_delta <- rep(1, nrow(train_s))
        }

        if(fluctuation == "logistic"){
          logitUpdate<- glm(selector[[v]][[s]]$train_s_nco[which(train_s$NCO_delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$NCObarAW[which(train_s$NCO_delta==1)])) +  selector[[v]][[s]]$H.AWnco[which(train_s$NCO_delta==1)], family='binomial', weights = selector[[v]][[s]]$wt_nco[which(train_s$NCO_delta==1)])

          epsilon <- logitUpdate$coef

          NCObarAW.star<- plogis(qlogis(selector[[v]][[s]]$NCObarAW) + epsilon*selector[[v]][[s]]$H.AWnco)
          NCObar1W.star<- plogis(qlogis(selector[[v]][[s]]$NCObar1W) + epsilon*selector[[v]][[s]]$H.1Wnco)
          NCObar0W.star<- plogis(qlogis(selector[[v]][[s]]$NCObar0W) + epsilon*selector[[v]][[s]]$H.0Wnco)

        } else {
          logitUpdate<- glm(selector[[v]][[s]]$train_s_nco[which(train_s$NCO_delta==1)] ~ -1 + offset(selector[[v]][[s]]$NCObarAW[which(train_s$NCO_delta==1)]) +  selector[[v]][[s]]$H.AWnco[which(train_s$NCO_delta==1)], family='gaussian', weights = selector[[v]][[s]]$wt_nco[which(train_s$NCO_delta==1)])

          epsilon <- logitUpdate$coef

          NCObarAW.star<- selector[[v]][[s]]$NCObarAW + epsilon*selector[[v]][[s]]$H.AWnco
          NCObar1W.star<- selector[[v]][[s]]$NCObar1W + epsilon*selector[[v]][[s]]$H.1Wnco
          NCObar0W.star<- selector[[v]][[s]]$NCObar0W + epsilon*selector[[v]][[s]]$H.0Wnco
        }

        if(family=="gaussian" & fluctuation == "logistic"){
          NCObarAW.star <- NCObarAW.star*(max(data$nco, na.rm = TRUE) - min(data$nco, na.rm = TRUE)) + min(data$nco, na.rm = TRUE)
          NCObar1W.star <- NCObar1W.star*(max(data$nco, na.rm = TRUE) - min(data$nco, na.rm = TRUE)) + min(data$nco, na.rm = TRUE)
          NCObar0W.star <- NCObar0W.star*(max(data$nco, na.rm = TRUE) - min(data$nco, na.rm = TRUE)) + min(data$nco, na.rm = TRUE)
          train_s_nco <- selector[[v]][[s]]$train_s_nco*(max(data$nco, na.rm = TRUE) - min(data$nco, na.rm = TRUE)) + min(data$nco, na.rm = TRUE)
        } else {
          train_s_nco <- selector[[v]][[s]]$train_s_nco
        }

        EICnco[which(data$v!=v & data$S %in% comparisons[[s]]),(length(comparisons)*(v-1)+s)] <- (selector[[v]][[s]]$wt_nco*selector[[v]][[s]]$H.AWnco*(train_s_nco - NCObarAW.star) + NCObar1W.star - NCObar0W.star - mean(NCObar1W.star - NCObar0W.star))/(length(which(train$S %in% comparisons[[s]]))/nrow(data))

        bias_nco[[v]][s] <- mean(NCObar1W.star - NCObar0W.star)

        addncobias[s] <- (bias[[v]][s] + bias_nco[[v]][s])^2 + var[s]
        ncoonly[s] <- (bias_nco[[v]][s])^2 + var[s]
      }
      b2v[s] <- bias[[v]][s]^2 + var[s]
    }

    lambdatilde$b2v[[v]] <- comparisons[[which(b2v==min(b2v))]]
    proportionselected$b2v[[v]] <- which(b2v==min(b2v))

    if(is.null(NCO)==FALSE){
      lambdatilde$ncobias[[v]] <- comparisons[[which(addncobias==min(addncobias))]]
      proportionselected$ncobias[[v]] <- which(addncobias==min(addncobias))

      lambdatilde$ncoonly[[v]] <- comparisons[[which(ncoonly==min(ncoonly))]]
      proportionselected$ncoonly[[v]] <- which(ncoonly==min(ncoonly))
    }

  }

  valid_initial <- validpreds(data, folds, V, selector, pRCT, Delta, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons)

  results$pooled$ATE <- list()
  results$pooled$ATE$b2v <- vector()
  results$pooled$ATE$ncobias <- vector()
  results$pooled$ATE$ncoonly <- vector()

  limitdist <- limitdistvar(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons)

  pool <- vector()
  for(v in 1:V){
    pool[v]<- limitdist$psi[[v]][proportionselected$b2v[[v]]]
  }
  results$pooled$ATE$b2v <- mean(pool)

  if(is.null(NCO)==FALSE){
    pool <- vector()
    for(v in 1:V){
      pool[v]<- limitdist$psi[[v]][proportionselected$ncobias[[v]]]
    }
    results$pooled$ATE$ncobias <- mean(pool)

    pool <- vector()
    for(v in 1:V){
      pool[v]<- limitdist$psi[[v]][proportionselected$ncoonly[[v]]]
    }
    results$pooled$ATE$ncoonly <- mean(pool)
  }

  #Estimated covariance matrices
  psipoundvec <- NA
  for(v in 1:V){
    psipoundvec <- c(psipoundvec,bias[[v]])
  }
  psipoundvec <- psipoundvec[-1]

  if(is.null(NCO)==FALSE){
    psipoundplusphivec <- NA
    for(v in 1:V){
      psipoundplusphivec <- c(psipoundplusphivec,(bias[[v]] + bias_nco[[v]]))
    }
    psipoundplusphivec <- psipoundplusphivec[-1]

    phivec <- NA
    for(v in 1:V){
      phivec <- c(phivec,bias_nco[[v]])
    }
    phivec <- phivec[-1]

    #make stacked ztilde vector

    #overall covariance matrix for ztilde_poundplusphi
    EICpoundplusphi <- EICpsipound+EICnco
    EICmat_poundplusphi <- cbind(EICpoundplusphi, limitdist$EICay)
    covMat_poundplusphi <- (t(EICmat_poundplusphi)%*%EICmat_poundplusphi)/nrow(data)

    #overall covariance matrix for ztilde_phi
    EICmat_phi <- cbind(EICnco, limitdist$EICay)
    covMat_phi <- (t(EICmat_phi)%*%EICmat_phi)/nrow(data)

    ztilde_poundplusphi_samp <- mvrnorm(n = 1000, mu=rep(0,ncol(EICmat_poundplusphi)), Sigma=covMat_poundplusphi/nrow(data))
    ztilde_phi_samp <- mvrnorm(n = 1000, mu=rep(0,ncol(EICmat_phi)), Sigma=covMat_phi/nrow(data))


  }

  #overall covariance matrix for ztilde
  EICmat <- cbind(EICpsipound, limitdist$EICay)

  covMat <- (t(EICmat)%*%EICmat)/nrow(data)

  #sample from multivariate ztildes
  ztilde_samp <- mvrnorm(n = 1000, mu=rep(0,ncol(EICmat)), Sigma=covMat/nrow(data))

  #selector for each sample
  biassample_psipound <- ztilde_samp[,(1:as.numeric(length(comparisons)*V))]
  if(is.null(NCO)==FALSE){
    biassample_psipoundplusphi <- ztilde_poundplusphi_samp[,(1:as.numeric(length(comparisons)*V))]
    biassample_phi <- ztilde_phi_samp[,(1:as.numeric(length(comparisons)*V))]
  }

  lambdatildeb2v <- matrix(NA, nrow=1000, ncol=length(psipoundvec))
  lambdatildencobias <- matrix(NA, nrow=1000, ncol=length(psipoundplusphivec))
  lambdatildencoonly <- matrix(NA, nrow=1000, ncol=length(phivec))
  for(b in 1:1000){
    lambdatildeb2v[b,] <- (biassample_psipound[b,]+psipoundvec)^2 + EICay
    if(is.null(NCO)==FALSE){
      lambdatildencobias[b,] <- (biassample_psipoundplusphi[b,] + psipoundplusphivec)^2 + EICay
      lambdatildencoonly[b,] <- (biassample_phi[b,] + phivec)^2 + EICay
    }
  }

  psisamp <- ztilde_samp[,(((as.numeric(length(comparisons)*V)+1)):(2*as.numeric(length(comparisons)*V)))]
  if(is.null(NCO)==FALSE){
    psisamp_poundplusphi <- ztilde_poundplusphi_samp[,(((as.numeric(length(comparisons)*V)+1)):(2*as.numeric(length(comparisons)*V)))]
    psisamp_phi <- ztilde_phi_samp[,(((as.numeric(length(comparisons)*V)+1)):(2*as.numeric(length(comparisons)*V)))]
  }

  #arrange V samples from limit distribution for psi_star for each sample
  sample_psi_pstarnv<- list()
  for(b in 1:1000){
    sample_psi_pstarnv[[b]] <- matrix(0, nrow=V, ncol=length(comparisons))
    for(v in 1:V){
      sample_psi_pstarnv[[b]][v,] <- psisamp[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]
    }
  }

  #now take average over whichever selected in the bias samples for each of 1000 samples
  psi_pstarnv_b2v <- vector()
  psi_pstarnv_b2v_v <- list()
  psi_pstarnv_nco <- vector()
  psi_pstarnv_nco_v <- list()
  psi_pstarnv_ncoonly <- vector()
  psi_pstarnv_ncoonly_v <- list()
  for(b in 1:1000){
    psi_pstarnv_b2v_v[[b]] <- vector()
    psi_pstarnv_nco_v[[b]] <- vector()
    psi_pstarnv_ncoonly_v[[b]] <- vector()
    for(v in 1:V){
      psi_pstarnv_b2v_v[[b]][v] <- sample_psi_pstarnv[[b]][v,which(lambdatildeb2v[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]==min(lambdatildeb2v[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]))]
      if(is.null(NCO)==FALSE){
        psi_pstarnv_nco_v[[b]][v] <- sample_psi_pstarnv[[b]][v,which(lambdatildencobias[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]==min(lambdatildencobias[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]))]
        psi_pstarnv_ncoonly_v[[b]][v] <- sample_psi_pstarnv[[b]][v,which(lambdatildencoonly[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]==min(lambdatildencoonly[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]))]
      }
    }
    psi_pstarnv_b2v[b] <- mean(psi_pstarnv_b2v_v[[b]])
    if(is.null(NCO)==FALSE){
      psi_pstarnv_nco[b] <- mean(psi_pstarnv_nco_v[[b]])
      psi_pstarnv_ncoonly[b] <- mean(psi_pstarnv_ncoonly_v[[b]])
    }
  }

  results$ATE <- list()
  results$Var <- list()

  results$ATE$b2v <- list()
  results$ATE$ncobias <- list()
  results$ATE$ncoonly <- list()

  results$Var$b2v <- list()
  results$Var$ncobias <- list()
  results$Var$ncoonly <- list()

  results$CI$b2v <- list()
  results$CI$ncobias <- list()
  results$CI$ncoonly <- list()

  results$ATE$b2v <- results$pooled$ATE$b2v
  if(any(unlist(proportionselected$b2v)!=1)){
    results$Var$b2v <- var(psi_pstarnv_b2v)
    results$CI$b2v <- results$ATE$b2v + quantile(psi_pstarnv_b2v, probs = c(0.025,0.5,0.975))

  } else {
    results$Var$b2v <- limitdist$Var
    results$CI$b2v <- c((results$pooled$ATE$b2v - 1.96*(limitdist$Var)^(1/2)), (results$pooled$ATE$b2v + 1.96*(limitdist$Var)^(1/2)))
  }


  if(is.null(NCO)==FALSE){
    results$ATE$ncobias <- results$pooled$ATE$ncobias
    if(any(unlist(proportionselected$ncobias)!=1)){
      results$Var$ncobias <- var(psi_pstarnv_nco)
      results$CI$ncobias <- results$ATE$ncobias + quantile(psi_pstarnv_nco, probs = c(0.025,0.5,0.975))

    } else {
      results$Var$ncobias <- limitdist$Var
      results$CI$ncobias <- c((results$pooled$ATE$ncobias - 1.96*(limitdist$Var)^(1/2)), (results$pooled$ATE$ncobias + 1.96*(limitdist$Var)^(1/2)))
    }

    results$ATE$ncoonly <- results$pooled$ATE$ncoonly
    if(any(unlist(proportionselected$ncoonly)!=1)){
      results$Var$ncoonly <- var(psi_pstarnv_ncoonly)
      results$CI$ncoonly <- results$ATE$ncoonly + quantile(psi_pstarnv_ncoonly, probs = c(0.025,0.5,0.975))

    } else {
      results$Var$ncoonly <- limitdist$Var
      results$CI$ncoonly <- c((results$pooled$ATE$ncoonly - 1.96*(limitdist$Var)^(1/2)), (results$pooled$ATE$ncoonly + 1.96*(limitdist$Var)^(1/2)))
    }
  }

  proportionselected_meanmat <- list()

  proportionselected_meanmat$b2v <- vector()
  if(is.null(NCO)==FALSE){
    proportionselected_meanmat$ncobias <- vector()
    proportionselected_meanmat$ncoonly<- vector()
  }

  results$lambdatilde <- list()
  results$lambdatilde$b2v <- vector()
  results$lambdatilde$ncobias <- vector()
  results$lambdatilde$ncoonly <- vector()

  for(v in 1:V){
    results$lambdatilde$b2v[v] <- max(lambdatilde$b2v[[v]])
    if(is.null(NCO)==FALSE){
      results$lambdatilde$ncobias[v] <- max(lambdatilde$ncobias[[v]])
      results$lambdatilde$ncoonly[v] <- max(lambdatilde$ncoonly[[v]])

    }
  }

  results$proportionselected_mean <- list()
  results$proportionselected_mean$b2v <- (mean(unlist(proportionselected$b2v))-1)
  if(is.null(NCO)==FALSE){
    results$proportionselected_mean$ncobias <- (mean(unlist(proportionselected$ncobias))-1)
    results$proportionselected_mean$ncoonly <- (mean(unlist(proportionselected$ncoonly))-1)
  }

  return(results)
}
