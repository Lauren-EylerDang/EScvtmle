#' @importFrom SuperLearner SuperLearner
#' @importFrom stats predict
#Function for training initial estimators for the outcome, treatment mechanism, and study selection mechanism regressions for RCT with or without real-world data (when active treatment is available in real-world data)
selector_func_txrwd <- function(train_s, data, Q.SL.library, d.SL.library.RCT, d.SL.library.RWD, g.SL.library, pRCT, family, family_nco, fluctuation = "logistic", NCO=NULL, Delta=NULL, Delta_NCO = NULL, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, bounds=NULL, cvControl=list()){

  #Train regressions for different experiments
  if(any(train_s$S==0)){

    Y <- train_s$Y
    if(family=="gaussian" & fluctuation == "logistic"){
      Y <- (Y - min(data$Y, na.rm = TRUE))/(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE))
    }

    if(adjustnco == FALSE){
      X <- train_s[,which((colnames(train_s) %in% c("Y", "nco", "NCO_delta", "v", "id"))==FALSE)]
    } else {
      X <- train_s[,which((colnames(train_s) %in% c("Y", "NCO_delta", "v", "id"))==FALSE)]
    }


    # set A=0 in X0 and A=1 in X1
    X0 <- X1 <- X
    X0$A <- 0
    X1$A <- 1

    Ynomiss <- Y[which(X$Delta==1)]
    Xnomiss <- X[which(X$Delta==1),]
    Xnomiss <- Xnomiss[ , -which(colnames(Xnomiss) %in% c("Delta"))]


    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL<- suppressWarnings(SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$Delta==1)]))
      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL$cvRisk)
        QbarSL <- QbarSL$fitLibrary[[keepAlg]]
      }

    } else {
      QbarSL<- suppressWarnings(SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, cvControl = cvControl, id=train_s$id[which(X$Delta==1)]))
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
      DbarSL<- suppressWarnings(SuperLearner(Y=X$Delta, X=X[ , -which(colnames(X) %in% c("Delta"))], SL.library=d.SL.library.RWD, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
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
    gHatSL<- suppressWarnings(SuperLearner(Y=X$A, X=X[ , -which(colnames(X) %in% c("A","Delta"))], SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
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

      # also want to evaluate the clever covariates at A=0 for all subjects
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

    #TMLE for E[E[Y|W,A=0,S=1]]

    # set S=1 and A to 0 or 1
    XS1 <- XS0 <- X
    XS1$S <- XS0$S <- 1
    XS1$A <- 1
    XS0$A <- 0

    # initial estimates of the outcome, given the observed exposure & covariates
    if(Q.discreteSL==TRUE){
      QbarSAW <- predict(QbarSL, newdata=X)
      # estimates of the outcome, given S=1 and A=1 or A=0 and covariates
      QbarS1W<- predict(QbarSL, newdata=XS1)
      QbarS0W<- predict(QbarSL, newdata=XS0)
    } else {
      QbarSAW <- predict(QbarSL, newdata=X)$pred
      # estimates of the outcome, given S=1 and A=1 or A=0 and covariates
      QbarS1W<- predict(QbarSL, newdata=XS1)$pred
      QbarS0W<- predict(QbarSL, newdata=XS0)$pred
    }

    # Estimate the trial participation mechanism g(S|A,W)
    gSHatSL<- suppressWarnings(SuperLearner(Y=X$S, X=X[ , -which(names(X) %in% c("S","A","Delta"))], SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
    if(g.discreteSL==TRUE){
      keepAlg <- which.min(gSHatSL$cvRisk)
      gSHatSL <- gSHatSL$fitLibrary[[keepAlg]]
      gSHat1W <- predict(gSHatSL, newdata = X)
    } else {
      # generate predicted prob of being in trial, given baseline covariates
      gSHat1W <- predict(gSHatSL, newdata = X)$pred
    }

    #get missingness mechanism setting S=1: P(Delta=1|A,S=1,W)
    if(is.null(Delta)==FALSE){
      if(d.discreteSL==TRUE){
        dHat1W$S1A1 <- predict(DbarSL, newdata = XS1)
        dHat1W$S1A0 <- predict(DbarSL, newdata = XS0)
      } else {
        dHat1W$S1A1 <- predict(DbarSL, newdata = XS1)$pred
        dHat1W$S1A0 <- predict(DbarSL, newdata = XS0)$pred
      }
    } else {
      dHat1W$S1A1 <- dHat1W$S1A0 <- rep(1, nrow(X))
    }

    #P(A|S=1,W)
    if(g.discreteSL==TRUE){
      gHata1_S1W<- predict(gHatSL, newdata = XS1)
    } else {
      gHata1_S1W<- predict(gHatSL, newdata = XS1)$pred
    }

    gHata0_S1W = 1-gHata1_S1W

    #-------------------------------------------------
    # Clever covariate H(S,A,W) for each subject
    #-------------------------------------------------
    if(target.gwt){
      wt_s <- as.numeric(X$A==1 & X$S==1 & X$Delta == 1)/.bound((gSHat1W*gHata1_S1W*dHat1W$S1A1),nrow(train_s), bounds) + as.numeric(X$A==0 & X$S==1 & X$Delta == 1)/.bound((gSHat1W*gHata0_S1W*dHat1W$S1A0),nrow(train_s), bounds)
      H.SAW1 <- as.numeric(X$A==1 & X$S==1 & X$Delta == 1)
      H.SAW0 <- -1*as.numeric(X$A==0 & X$S==1 & X$Delta == 1)

      # also want to evaluate the clever covariates at S=1 and A=1 or 0 for all subjects
      H.S1W<- rep(1, nrow(X))
      H.S0W<- rep(-1, nrow(X))
    } else{
      wt_s <- rep(1, nrow(X))
      H.SAW1 <- as.numeric(X$A==1 & X$S==1 & X$Delta == 1)/.bound((gSHat1W*gHata1_S1W*dHat1W$S1A1),nrow(train_s), bounds)
      H.SAW0 <- -1*as.numeric(X$A==0 & X$S==1 & X$Delta == 1)/.bound((gSHat1W*gHata0_S1W*dHat1W$S1A0),nrow(train_s), bounds)

      # also want to evaluate the clever covariates at S=1 and A=1 or 0 for all subjects
      H.S1W<- 1/.bound((gSHat1W*gHata1_S1W*dHat1W$S1A1),nrow(train_s), bounds)
      H.S0W<- -1/.bound((gSHat1W*gHata0_S1W*dHat1W$S1A0),nrow(train_s), bounds)
    }

    if(is.null(NCO)==FALSE){
      # call Super Learner for estimation of NCObarAW
      if(family_nco=="gaussian" & fluctuation == "logistic"){
        train_s_nco <- (train_s$nco - min(data$nco, na.rm=TRUE))/(max(data$nco, na.rm=TRUE) - min(data$nco, na.rm=TRUE))
      } else {
        train_s_nco <- train_s$nco
      }

      X <- train_s[,which((colnames(train_s) %in% c("Y", "nco", "Delta", "v", "id"))==FALSE)]
      X0 <- X1 <- X
      X0$A <- 0
      X1$A <- 1


      NCOnomiss <- train_s_nco[which(X$NCO_delta==1)]
      Xnomiss <- X[which(X$NCO_delta==1),]
      Xnomiss <- Xnomiss[ , -which(colnames(Xnomiss) %in% c("NCO_delta"))]


      if(fluctuation == "logistic"){
        NCObarSL<- suppressWarnings(SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$NCO_delta==1)]))
        if(Q.discreteSL==TRUE){
          keepAlg <- which.min(NCObarSL$cvRisk)
          NCObarSL <- NCObarSL$fitLibrary[[keepAlg]]
        }
      } else {
        NCObarSL<- suppressWarnings(SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$NCO_delta==1)]))
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
      gHatSLnco<- suppressWarnings(SuperLearner(Y=X$A, X=X[ , -which(colnames(X) %in% c("A","NCO_delta"))], SL.library=g.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
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
        DbarSLnco<- suppressWarnings(SuperLearner(Y=X$NCO_delta, X=X[ , -which(colnames(X) %in% c("NCO_delta"))], SL.library=d.SL.library.RWD, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
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
      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "nco", "NCO_delta", "v", "id"))==FALSE)]
    } else {
      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "NCO_delta", "v", "id"))==FALSE)]
    }


    # set the A=0 in X0, A=1 in X1
    X0 <- X1 <- X
    X0$A <- 0
    X1$A <- 1

    Ynomiss <- Y[which(X$Delta==1)]
    Xnomiss <- X[which(X$Delta==1),]
    Xnomiss <- Xnomiss[ , -which(colnames(Xnomiss) %in% c("Delta"))]

    # call Super Learner for estimation of QbarAW
    if(fluctuation == "logistic"){
      QbarSL<- suppressWarnings(SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$Delta==1)]))
      if(Q.discreteSL==TRUE){
        keepAlg <- which.min(QbarSL$cvRisk)
        QbarSL <- QbarSL$fitLibrary[[keepAlg]]
      }
    } else {
      QbarSL<- suppressWarnings(SuperLearner(Y=Ynomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$Delta==1)]))
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
      DbarSL<- suppressWarnings(SuperLearner(Y=X$Delta, X=X[ , -which(colnames(X) %in% c("Delta"))], SL.library=d.SL.library.RCT, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
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

    # Known probability of being exposed for RCT
    gHatSL<- NULL
    gHat1W<- rep(pRCT, nrow(X))
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

      X <- train_s[,which((colnames(train_s) %in% c("S","Y", "nco", "Delta", "v", "id"))==FALSE)]
      X0 <- X1 <- X
      X0$A <- 0
      X1$A <- 1

      NCOnomiss <- train_s_nco[which(X$NCO_delta==1)]
      Xnomiss <- X[which(X$NCO_delta==1),]
      Xnomiss <- Xnomiss[ , -which(colnames(Xnomiss) %in% c("NCO_delta"))]


      if(fluctuation == "logistic"){
        NCObarSL<- suppressWarnings(SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$NCO_delta==1)]))
        if(Q.discreteSL==TRUE){
          keepAlg <- which.min(NCObarSL$cvRisk)
          NCObarSL <- NCObarSL$fitLibrary[[keepAlg]]
        }
      } else {
        NCObarSL<- suppressWarnings(SuperLearner(Y=NCOnomiss, X=Xnomiss, SL.library=Q.SL.library, family=family, control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id[which(X$NCO_delta==1)]))
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
        DbarSLnco<- suppressWarnings(SuperLearner(Y=X$NCO_delta, X=X[ , -which(colnames(X) %in% c("NCO_delta"))], SL.library=d.SL.library.RCT, family="binomial", control = list(saveFitLibrary=TRUE), cvControl = cvControl, id=train_s$id))
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

        # also want to evaluate the clever covariates at A=0 and A=1 for all subjects
        H.0Wnco<- rep(-1, nrow(X))
        H.1Wnco<- rep(1, nrow(X))
      } else{
        wt_nco <- rep(1, nrow(X))
        H.AW1nco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)
        H.AW0nco <- -1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)

        H.AWnco <- (as.numeric(X$A==1 & X$NCO_delta==1))/.bound((gHat1W*dHat1Wnco$A1),nrow(train_s), bounds)-1*(as.numeric(X$A==0 & X$NCO_delta==1))/.bound((gHat0W*dHat1Wnco$A0),nrow(train_s), bounds)

        # also want to evaluate the clever covariates at A=0 and A=1 for all subjects
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

  out <- list("Y" = Y, "QbarSL" = QbarSL, "QbarAW" = QbarAW, "Qbar1W" = Qbar1W, "Qbar0W" = Qbar0W, "QbarSAW" = QbarSAW, "QbarS1W" = QbarS1W, "QbarS0W" = QbarS0W, "gHatSL" = gHatSL, "train_s_nco" = train_s_nco, "NCObarAW" = NCObarAW, "NCObar1W" = NCObar1W, "NCObar0W" = NCObar0W, "H.AW1" = H.AW1, "H.AW0" = H.AW0,"H.AW" = H.AW, "H.0W" = H.0W, "H.1W" = H.1W, "H.SAW1" = H.SAW1,"H.SAW0" = H.SAW0, "H.S1W" = H.S1W, "H.S0W" = H.S0W, "DbarSL" = DbarSL, "H.AWnco" = H.AWnco, "H.0Wnco" = H.0Wnco, "H.1Wnco"= H.1Wnco, "wt"=wt, "wt_s"=wt_s, "wt_nco"=wt_nco)

  return(out)

}

#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom stats qlogis
#Function to estimate the bias-variance tradeoff using the experiment-selection set for each fold when active treatment is available in the real-world data
bvt_txinrwd <- function(v, selector, NCO, comparisons, train, data, fluctuation, family, n.id){
  out <- list()

  out$b2v <- vector()
  out$bias <- vector()
  out$var <- vector()
  out$EIClambdav <- list()

  datid <- data[which(duplicated(data$id)==FALSE),]

  if(is.null(NCO)==FALSE){
    out$addncobias <- vector()
    out$bias_nco <- vector()
  }

  out$EICpsipound <- matrix(0, nrow=n.id, ncol=length(comparisons))
  out$EICnco <- matrix(0, nrow=n.id, ncol=length(comparisons))


  for(s in 1:length(comparisons)){
    train_s <- train[which(train$S %in% comparisons[[s]]),]

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

    #EIC
    EIClambdav_s <- (selector[[v]][[s]]$wt*(selector[[v]][[s]]$H.AW1 + selector[[v]][[s]]$H.AW0)*(trainsY - QbarAW.star) + Qbar1W.star - Qbar0W.star - mean(Qbar1W.star - Qbar0W.star))
    #account for dependent units
    #take mean of EIC observations by trains_s$id (confirmed taking mean correctly)
    out$EIClambdav[[s]] <- as.vector(by(EIClambdav_s, train_s$id, mean))

    if(s==1){
      QbarS1W.star <- Qbar1W.star
      QbarS0W.star <- Qbar0W.star
    } else {
      # Update the initial estimator
      if(fluctuation == "logistic"){
        logitUpdate<- glm(selector[[v]][[s]]$Y[which(train_s$Delta==1)] ~ -1 + offset(qlogis(selector[[v]][[s]]$QbarSAW[which(train_s$Delta==1)])) +  selector[[v]][[s]]$H.SAW1[which(train_s$Delta==1)] + selector[[v]][[s]]$H.SAW0[which(train_s$Delta==1)], family='quasibinomial', weights = selector[[v]][[s]]$wt_s[which(train_s$Delta==1)])

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

      if(family=="gaussian" & fluctuation == "logistic"){
        QbarSAW.star <- QbarSAW.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        QbarS1W.star <- QbarS1W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
        QbarS0W.star <- QbarS0W.star*(max(data$Y, na.rm = TRUE) - min(data$Y, na.rm = TRUE)) + min(data$Y, na.rm = TRUE)
      }
    }

    if(s>1){
      EICpsipound_s <- (selector[[v]][[s]]$wt*(selector[[v]][[s]]$H.AW1 + selector[[v]][[s]]$H.AW0)*(trainsY - QbarAW.star) + Qbar1W.star - Qbar0W.star - mean(Qbar1W.star - Qbar0W.star) - (selector[[v]][[s]]$wt_s*(selector[[v]][[s]]$H.SAW1 + selector[[v]][[s]]$H.SAW0)*(trainsY - QbarSAW.star) + QbarS1W.star - QbarS0W.star - mean(QbarS1W.star - QbarS0W.star)))
      EICpsipound_s <- as.vector(by(EICpsipound_s, train_s$id, mean))
      out$EICpsipound[which(datid$v!=v & datid$S %in% comparisons[[s]]),s] <- EICpsipound_s/(length(which(train$S %in% comparisons[[s]] & duplicated(train$id)==FALSE))/n.id)
    }

    out$bias[s] <- mean(Qbar1W.star - Qbar0W.star) - mean(QbarS1W.star - QbarS0W.star)
    out$var[s] <- var(out$EIClambdav[[s]])/length(unique(train_s$id))

    #nco
    if(is.null(NCO)==FALSE){

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

      EICnco_s <- as.vector(selector[[v]][[s]]$wt_nco*selector[[v]][[s]]$H.AWnco*(train_s_nco - NCObarAW.star) + NCObar1W.star - NCObar0W.star - mean(NCObar1W.star - NCObar0W.star))
      EICnco_s <- as.vector(by(EICnco_s, train_s$id, mean))

      out$EICnco[which(datid$v!=v & datid$S %in% comparisons[[s]]),s] <- EICnco_s/(length(which(train$S %in% comparisons[[s]] & duplicated(train$id)==FALSE))/n.id)

      out$bias_nco[s] <- mean(NCObar1W.star - NCObar0W.star)

      out$addncobias[s] <- (out$bias[s] + out$bias_nco[s])^2 + out$var[s]
    }
    out$b2v[s] <- out$bias[s]^2 + out$var[s]
  }
  return(out)
}
