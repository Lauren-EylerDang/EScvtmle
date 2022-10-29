#' @importFrom SuperLearner SuperLearner
#' @importFrom stats predict
#Function for training initial estimators for the outcome, treatment mechanism, and study selection mechanism regressions for RCT with or without real-world data (when no active treatment is available in real-world data)
selector_func_notxrwd <- function(train_s, data, Q.SL.library, d.SL.library, g.SL.library, pRCT = pRCT, family, family_nco, fluctuation = "logistic", NCO=NULL, Delta=NULL, Delta_NCO = NULL, adjustnco=adjustnco, target.gwt=target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL, bounds){

  #Train regressions for different experiments
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


    # set A=0 in X0 and A=1 in X1
    X0 <- X1 <- X
    X0$A <- 0
    X1$A <- 1

    if(is.null(Delta)==FALSE){
      Ynomiss <- Y[which(X$Delta==1)]
      Xnomiss <- X[which(X$Delta==1),]
      Xnomiss <- Xnomiss[ , -which(colnames(Xnomiss) %in% c("Delta"))]
    } else {
      Ynomiss <- Y
      Xnomiss <- X
    }

    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL<- suppressWarnings(SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE)))
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
      # estimates of the outcome, given A=0 or A=1 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)
      Qbar1W<- predict(QbarSL, newdata=X1)
    } else {
      QbarAW <- predict(QbarSL, newdata=X)$pred
      # estimates of the outcome, given A=0 or A=1 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)$pred
      Qbar1W<- predict(QbarSL, newdata=X1)$pred
    }

    dHat1W <- list()
    if(is.null(Delta)==FALSE){
      DbarSL<- SuperLearner(Y=X$Delta, X=X[ , -which(colnames(X) %in% c("Delta"))], SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
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
    gHatSL<- SuperLearner(Y=X$A, X=X[ , -which(colnames(X) %in% c("A","Delta"))], SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
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
      wt <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds) + 1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)
      H.AW1 <- as.numeric(X$A==1 & X$Delta==1)
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))
      H.AW <- (as.numeric(X$A==1 & X$Delta==1))-1*(as.numeric(X$A==0 & X$Delta==1))

      # also want to evaluate the clever covariates at A=0 or A=1 for all subjects
      H.0W<- rep(-1, nrow(X))
      H.1W<- rep(1, nrow(X))
    } else{
      wt <- rep(1, nrow(X))
      H.AW1 <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds)
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)

      H.AW <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds)-1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)

      # also want to evaluate the clever covariates at A=0 or A=1 for all subjects
      H.0W<- -1/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)
      H.1W<- 1/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds)
    }

    #TMLE for E[E[Y|W,A=0,S=1]]
    if(adjustnco == FALSE){
      XS <- train_s[,which((colnames(train_s) %in% c("Y", "nco", "NCO_delta", "v"))==FALSE)]
    } else {
      XS <- train_s[,which((colnames(train_s) %in% c("Y", "NCO_delta", "v"))==FALSE)]
    }


    # set A=0, S=1
    XS0 <-XS
    XS0$A <- 0 # under control
    XS0$S <- 1

    if(is.null(Delta)==FALSE){
      XSnomiss <- XS[which(XS$Delta==1),]
      XSnomiss <- XSnomiss[ , -which(colnames(XSnomiss) %in% c("Delta"))]
    } else {
      XSnomiss <- XS
      XS$Delta <- rep(1, nrow(XS))
    }
    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL_S<- suppressWarnings(SuperLearner(Y=Ynomiss, X=XSnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE)))

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
      # estimates of the outcome, given A=0, S=1, and covariates
      QbarS0W<- predict(QbarSL_S, newdata=XS0)
    } else {
      QbarSAW <- predict(QbarSL_S, newdata=XS)$pred
      # estimates of the outcome, given A=0, S=1, and covariates
      QbarS0W<- predict(QbarSL_S, newdata=XS0)$pred
    }

    # Estimate the trial participation mechanism g(S|A=0,Delta=1,W)
    #------------------------------------------
    lambda_A0 <- XS[which(XS$A==0 & XS$Delta==1),]
    gSHatSL<- SuperLearner(Y=lambda_A0$S, X=lambda_A0[ , -which(colnames(lambda_A0) %in% c("A","S","Delta"))], SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
    if(g.discreteSL==TRUE){
      keepAlg <- which.min(gSHatSL$cvRisk)
      gSHatSL <- gSHatSL$fitLibrary[[keepAlg]]
      gSHat1W <- predict(gSHatSL, newdata = XS)
    } else {
      gSHat1W <- predict(gSHatSL, newdata = XS)$pred
    }


    #-------------------------------------------------
    # Clever covariate H(S,A,W) for each subject
    #-------------------------------------------------
    if(target.gwt){
      wt_s <- as.numeric(XS$A==0 & XS$S==1 & XS$Delta == 1)/.bound((gSHat1W*gHat0W*dHat1W$A0),nrow(train_s), bounds)
      H.SAW<- as.numeric(XS$A==0 & XS$S==1 & XS$Delta == 1)

      # also want to evaluate the clever covariates at S=1 and A=0 for all subjects
      H.S0W<- rep(1, nrow(XS))
    } else{
      wt_s <- rep(1, nrow(XS))
      H.SAW<- as.numeric(XS$A==0 & XS$S==1 & XS$Delta == 1)/.bound((gSHat1W*gHat0W*dHat1W$A0),nrow(train_s), bounds)

      # also want to evaluate the clever covariates at S=1 and A=0 for all subjects
      H.S0W<- 1/.bound((gSHat1W*gHat0W*dHat1W$A0),nrow(train_s), bounds)
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
        Xnomiss <- Xnomiss[ , -which(names(Xnomiss) %in% c("NCO_delta"))]
      } else {
        NCOnomiss <- train_s_nco
        Xnomiss <- X
        X$NCO_delta <- rep(1, nrow(X))
      }

      if(fluctuation == "logistic"){
        NCObarSL<- suppressWarnings(SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE)))
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
      gHatSLnco<- SuperLearner(Y=X$A, X=X[ , -which(colnames(X) %in% c("A","NCO_delta"))], SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
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
        DbarSLnco<- SuperLearner(Y=X$NCO_delta, X=X[ , -which(colnames(X) %in% c("NCO_delta"))], SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
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
        wt_nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds) + 1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)
        H.AW1nco <- as.numeric(X$A==1 & X$NCO_delta==1)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))
        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))-1*(as.numeric(X$A==0 & X$NCO_delta==1))

        # also want to evaluate the clever covariates at A=0 for all subjects
        H.0Wnco<- rep(-1, nrow(X))
        H.1Wnco<- rep(1, nrow(X))
      } else{
        wt_nco <- rep(1, nrow(X))
        H.AW1nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)

        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)-1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)

        # also want to evaluate the clever covariates at A=0 for all subjects
        H.0Wnco<- -1/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)
        H.1Wnco<- 1/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)
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


    # set A=0 in X0 and A=1 in X1
    X0 <- X1 <- X
    X0$A <- 0
    X1$A <- 1

    if(is.null(Delta)==FALSE){
      Ynomiss <- Y[which(X$Delta==1)]
      Xnomiss <- X[which(X$Delta==1),]
      Xnomiss <- Xnomiss[ , -which(colnames(Xnomiss) %in% c("Delta"))]
    } else {
      Ynomiss <- Y
      Xnomiss <- X
    }

    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL<- suppressWarnings(SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE)))
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
      # estimates of the outcome, given A=0 or A=1 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)
      Qbar1W<- predict(QbarSL, newdata=X1)
    } else {
      # initial estimates of the outcome, given the observed exposure & covariates
      QbarAW <- predict(QbarSL, newdata=X)$pred
      # estimates of the outcome, given A=0 or A=1 and covariates
      Qbar0W<- predict(QbarSL, newdata=X0)$pred
      Qbar1W<- predict(QbarSL, newdata=X1)$pred
    }

    dHat1W <- list()
    if(is.null(Delta)==FALSE){
      DbarSL<- SuperLearner(Y=X$Delta, X=X[ , -which(colnames(X) %in% c("Delta"))], SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
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

    # Known randomization probability for RCT
    gHatSL<- NULL
    gHat1W<- rep(pRCT, nrow(X))
    gHat0W<- 1- gHat1W

    #-------------------------------------------------
    # Clever covariate H(A,W) for each subject
    #-------------------------------------------------
    if(target.gwt){
      wt <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds) + 1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)
      H.AW1 <- as.numeric(X$A==1 & X$Delta==1)
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))
      H.AW <- (as.numeric(X$A==1 & X$Delta==1))-1*(as.numeric(X$A==0 & X$Delta==1))

      # also want to evaluate the clever covariates at A=0 and A=1 for all subjects
      H.0W<- rep(-1, nrow(X))
      H.1W<- rep(1, nrow(X))
    } else{
      wt <- rep(1, nrow(X))
      H.AW1 <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds)
      H.AW0 <- -1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)

      H.AW <- (as.numeric(X$A==1 & X$Delta==1))/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds)-1*(as.numeric(X$A==0 & X$Delta==1))/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)

      # also want to evaluate the clever covariates at A=0 for all subjects
      H.0W<- -1/.bound((gHat0W*dHat1W$A0),nrow(train_s), bounds)
      H.1W<- 1/.bound((gHat1W*dHat1W$A1),nrow(train_s), bounds)
    }

    QbarSAW <- NULL
    QbarS0W <- NULL
    H.SAW <- NULL
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
        Xnomiss <- Xnomiss[ , -which(names(Xnomiss) %in% c("NCO_delta"))]
      } else {
        NCOnomiss <- train_s_nco
        Xnomiss <- X
      }

      if(fluctuation == "logistic"){
        NCObarSL<- suppressWarnings(SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE)))
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
      gHat0W<- 1- gHat1W

      dHat1Wnco <- list()
      if(is.null(Delta_NCO)==FALSE){
        DbarSLnco<- SuperLearner(Y=X$NCO_delta, X=X[ , -which(colnames(X) %in% c("NCO_delta"))], SL.library=d.SL.library, family="binomial", control = list(saveFitLibrary=TRUE))
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
        wt_nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds) + 1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)
        H.AW1nco <- as.numeric(X$A==1 & X$NCO_delta==1)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))
        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))-1*(as.numeric(X$A==0 & X$NCO_delta==1))

        # also want to evaluate the clever covariates at A=0 or A=1 for all subjects
        H.0Wnco<- rep(-1, nrow(X))
        H.1Wnco<- rep(1, nrow(X))
      } else{
        wt_nco <- rep(1, nrow(X))
        H.AW1nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)

        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)-1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)

        # also want to evaluate the clever covariates at A=0 or A=1 for all subjects
        H.0Wnco<- -1/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)
        H.1Wnco<- 1/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)
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

  out <- list("Y" = Y, "QbarSL" = QbarSL, "QbarAW" = QbarAW, "Qbar1W" = Qbar1W, "Qbar0W" = Qbar0W, "QbarSAW" = QbarSAW, "QbarS0W" = QbarS0W, "gHatSL" = gHatSL, "train_s_nco" = train_s_nco, "NCObarAW" = NCObarAW, "NCObar1W" = NCObar1W, "NCObar0W" = NCObar0W, "H.AW1" = H.AW1, "H.AW0" = H.AW0,"H.AW" = H.AW, "H.0W" = H.0W, "H.1W" = H.1W, "H.SAW" = H.SAW, "H.S0W" = H.S0W, "DbarSL" = DbarSL, "H.AWnco" = H.AWnco, "H.0Wnco" = H.0Wnco, "H.1Wnco"= H.1Wnco, "wt"=wt, "wt_s"=wt_s, "wt_nco"=wt_nco)

  return(out)

}


#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom stats qlogis
#Function to estimate the bias-variance tradeoff using the experiment-selection set for each fold when active treatment is not available in the real-world data
bvt_notxinrwd <- function(v, selector, NCO, comparisons, train, data, fluctuation, family){
  out <- list()

  out$b2v <- vector()
  out$bias <- vector()
  out$var <- vector()
  out$EIClambdav <- list()

  if(is.null(NCO)==FALSE){
    out$addncobias <- vector()
    out$bias_nco <- vector()
  }

  out$EICpsipound <- matrix(0, nrow=nrow(data), ncol=length(comparisons))
  out$EICnco <- matrix(0, nrow=nrow(data), ncol=length(comparisons))


  for(s in 1:length(comparisons)){
    train_s <- train[which(train$S %in% comparisons[[s]]),]
    if(is.null(train_s$Delta)){
      train_s$Delta <- rep(1, nrow(train_s))
    }

    #tmle
    if(fluctuation == "logistic"){
      logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$QbarAW[which(train_s$Delta==1)])) +  selector[[v]][[s]]$H.AW1[which(train_s$Delta==1)] + selector[[v]][[s]]$H.AW0[which(train_s$Delta==1)], family='quasibinomial', weights = selector[[v]][[s]]$wt[which(train_s$Delta==1)])

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

    out$EIClambdav[[s]] <- (selector[[v]][[s]]$wt*(selector[[v]][[s]]$H.AW1 + selector[[v]][[s]]$H.AW0)*(trainsY - QbarAW.star) + Qbar1W.star - Qbar0W.star - mean(Qbar1W.star - Qbar0W.star))

    if(s==1){
      QbarS0W.star <- Qbar0W.star
    } else {
      # Update the initial estimator
      if(fluctuation == "logistic"){
        logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$QbarSAW[which(train_s$Delta==1)])) +  selector[[v]][[s]]$H.SAW[which(train_s$Delta==1)], family='quasibinomial', weights = selector[[v]][[s]]$wt_s[which(train_s$Delta==1)])

        epsilon <- logitUpdate$coef

        QbarSAW.star<- plogis(qlogis(selector[[v]][[s]]$QbarSAW) + epsilon*selector[[v]][[s]]$H.SAW)
        QbarS0W.star<- plogis(qlogis(selector[[v]][[s]]$QbarS0W) + epsilon*selector[[v]][[s]]$H.S0W)
      } else {
        logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(selector[[v]][[s]]$QbarSAW[which(train_s$Delta==1)]) +  selector[[v]][[s]]$H.SAW[which(train_s$Delta==1)], family='gaussian', weights = selector[[v]][[s]]$wt_s[which(train_s$Delta==1)])

        epsilon <- logitUpdate$coef

        QbarSAW.star<- selector[[v]][[s]]$QbarSAW + epsilon*selector[[v]][[s]]$H.SAW
        QbarS0W.star<- selector[[v]][[s]]$QbarS0W + epsilon*selector[[v]][[s]]$H.S0W
      }

      if(family=="gaussian" & fluctuation == "logistic"){
        QbarSAW.star <- QbarSAW.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        QbarS0W.star <- QbarS0W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
      }
    }

    if(s>1){
      out$EICpsipound[which(data$v!=v & data$S %in% comparisons[[s]]),s] <- (selector[[v]][[s]]$wt_s*(selector[[v]][[s]]$H.SAW)*(trainsY - QbarSAW.star) + QbarS0W.star - mean(QbarS0W.star) + selector[[v]][[s]]$wt*(selector[[v]][[s]]$H.AW0)*(trainsY - QbarAW.star) - Qbar0W.star + mean(Qbar0W.star))/(length(which(train$S %in% comparisons[[s]]))/nrow(data))
    }

    out$bias[s] <- mean(QbarS0W.star) - mean(Qbar0W.star)
    out$var[s] <- var(out$EIClambdav[[s]])/length(trainsY)

    #nco
    if(is.null(NCO)==FALSE){

      if(is.null(train_s$NCO_delta)){
        train_s$NCO_delta <- rep(1, nrow(train_s))
      }

      if(fluctuation == "logistic"){
        logitUpdate<- glm(selector[[v]][[s]]$train_s_nco[which(train_s$NCO_delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$NCObarAW[which(train_s$NCO_delta==1)])) +  selector[[v]][[s]]$H.AWnco[which(train_s$NCO_delta==1)], family='quasibinomial', weights = selector[[v]][[s]]$wt_nco[which(train_s$NCO_delta==1)])

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

      out$EICnco[which(data$v!=v & data$S %in% comparisons[[s]]),s] <- (selector[[v]][[s]]$wt_nco*selector[[v]][[s]]$H.AWnco*(train_s_nco - NCObarAW.star) + NCObar1W.star - NCObar0W.star - mean(NCObar1W.star - NCObar0W.star))/(length(which(train$S %in% comparisons[[s]]))/nrow(data))

      out$bias_nco[s] <- mean(NCObar1W.star - NCObar0W.star)

      out$addncobias[s] <- (out$bias[s] + out$bias_nco[s])^2 + out$var[s]
    }
    out$b2v[s] <- out$bias[s]^2 + out$var[s]
  }
  return(out)
}

