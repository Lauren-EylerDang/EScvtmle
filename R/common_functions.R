#'
#bound denominator of clever covariates
.bound <- function(x,n, bounds=NULL){
  if(is.null(bounds)){
    x <- pmax((5/(n)^(1/2)/log(n)), pmin(1,x))
  } else {
    x <- pmax(bounds[1], pmin(bounds[2],x))
  }
  return(x)
}

#' @importFrom dplyr rename
#' @importFrom tidyselect all_of
# Function to pre-process data
# Includes removal of observations from observational dataset with W covariates not represented in RCT if txinrwd==FALSE
preprocess <- function(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, adjustnco = adjustnco){

  #remove observations missing treatment
  data <- rename(data, A = all_of(treatment_var))

  if(length(which(is.na(data$A)))>0) message("Removing observations with missing treatment variable.")
  data <- data[which(is.na(data$A)==FALSE),]

  #make A coded as 1=treatment of interest, 0=control
  data$A <- ifelse(data$A == treatment, 1, 0)

  data <- rename(data, S = all_of(study))
  data <- rename(data, Y = all_of(outcome))

  if (length(which(data$S>1 & data$A==1))>0 & txinrwd==FALSE) stop("Active treatment available in external data. Set txinrwd==TRUE.")

  #trim data to avoid positivity violation if txinrwd=FALSE
  if(txinrwd==FALSE){
    for(w in 1:length(covariates)){
      if(is.factor(data[,covariates[w]])==FALSE){
        whichW <- which(data$S!=1 & (data[,covariates[w]] < (min(data[which(data$S==1),covariates[w]]))) | (data[,covariates[w]] > (max(data[which(data$S==1),covariates[w]]))))
      } else {
        whichW <- which(data$S!=1 & (data[,covariates[w]] %in% (data[which(data$S==1),covariates[w]]))==FALSE)
      }
      if(length(whichW)>0){
        data <- data[-whichW,]
      }
    }
  }

  if(is.null(NCO) == FALSE){
    data <- rename(data, nco = all_of(NCO))

    if(adjustnco == TRUE & txinrwd == FALSE){
      whichW <- which(data$S!=1 & (data[,"nco"] < (min(data[which(data$S==1),"nco"]))) | (data[,"nco"] > (max(data[which(data$S==1),"nco"]))))
      if(length(whichW)>0){
        data <- data[-whichW,]
      }
    }

    if(is.null(Delta_NCO)==FALSE){
      data <- rename(data, NCO_delta = all_of(Delta_NCO))
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
    data <- rename(data, Delta = all_of(Delta))
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

  return(data)
}


#apply selector_func to different datasets
apply_selector_func <- function(txinrwd, train, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO=NULL, Delta=NULL, Delta_NCO=NULL, adjustnco=adjustnco, target.gwt=target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL, comparisons, bounds){
  out <- list()
  for(s in 1:(length(comparisons))){

    train_s <- train[which(train$S %in% comparisons[[s]]),]
    train_s$S[which(train_s$S!=1)]<-0

    if(txinrwd==TRUE){
      out[[s]] <- selector_func_txrwd(train_s, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, bounds)
    } else {
      out[[s]] <- selector_func_notxrwd(train_s, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, bounds)
    }
  }
  return(out)
}

#' @importFrom stats predict
# Get initial estimates of the conditional mean outcome and treatment mechanism for validation set observations using regressions trained on training sets
validpreds <- function(data, folds, V, selector, pRCT, Delta, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons){
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
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$QbarAW <- predict(selector[[v]][[s]]$QbarSL, newdata = data[which(data$v==v & (data$S %in% comparisons[[s]])),])
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar1W <- predict(selector[[v]][[s]]$QbarSL, newdata = D1[which(data$v==v & (data$S %in% comparisons[[s]])),])
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar0W <- predict(selector[[v]][[s]]$QbarSL, newdata = D0[which(data$v==v & (data$S %in% comparisons[[s]])),])
      } else {
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$QbarAW <- predict(selector[[v]][[s]]$QbarSL, newdata = data[which(data$v==v & (data$S %in% comparisons[[s]])),])$pred
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar1W <- predict(selector[[v]][[s]]$QbarSL, newdata = D1[which(data$v==v & (data$S %in% comparisons[[s]])),])$pred
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$Qbar0W <- predict(selector[[v]][[s]]$QbarSL, newdata = D0[which(data$v==v & (data$S %in% comparisons[[s]])),])$pred
      }

      if(is.null(Delta)==FALSE){
        if(d.discreteSL==TRUE){
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbarAW <- predict(selector[[v]][[s]]$DbarSL, newdata = data[which(data$v==v & (data$S %in% comparisons[[s]])),])
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar1W <- predict(selector[[v]][[s]]$DbarSL, newdata = D1[which(data$v==v & (data$S %in% comparisons[[s]])),])
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar0W <- predict(selector[[v]][[s]]$DbarSL, newdata = D0[which(data$v==v & (data$S %in% comparisons[[s]])),])
        } else {
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbarAW <- predict(selector[[v]][[s]]$DbarSL, newdata = data[which(data$v==v & (data$S %in% comparisons[[s]])),])$pred
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar1W <- predict(selector[[v]][[s]]$DbarSL, newdata = D1[which(data$v==v & (data$S %in% comparisons[[s]])),])$pred
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$dbar0W <- predict(selector[[v]][[s]]$DbarSL, newdata = D0[which(data$v==v & (data$S %in% comparisons[[s]])),])$pred
        }
      }

      if(s==1){
        out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$gHat1W <- rep(pRCT, length(which(data$v==v & (data$S %in% comparisons[[s]]))))
      } else {
        if(length(comparisons)>1){
          out[[s]][which(data$v==v & (data$S %in% comparisons[[s]])),]$gHat1W <- predict(selector[[v]][[s]]$gHatSL, newdata = data[which(data$v==v & (data$S %in% comparisons[[s]])),])
        }
      }
    }
  }
  return(out)
}

#' @importFrom stringr str_match
#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom stats qlogis
#function to estimate components of limit distribution of ES-CVTMLE estimator
limitdistvar<- function(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons, bounds){
  out <- list()

  out$EICay <- matrix(0, nrow=nrow(data), ncol=length(comparisons)*V)
  out$psi <- list()
  for(v in 1:V){
    out$psi[[v]] <- vector()
  }

  out$clevercov <- list()

  for(s in 1:length(comparisons)){

    validmat <- valid_initial[[s]]
    validmat$Yscale <- data$Y

    if(family=="gaussian" & fluctuation == "logistic"){
      validmat$Yscale <- (validmat$Yscale - min(data$Y))/(max(data$Y) - min(data$Y))
    }

    if(("Delta" %in% colnames(data))==FALSE){
      data$Delta <- rep(1, nrow(data))
    }

    if(target.gwt){
      wt <- as.numeric(data$A==1 & data$Delta==1)/.bound((validmat$gHat1W*validmat$dbar1W),length(which(data$S %in% comparisons[[s]])), bounds) + as.numeric(data$A==0 & data$Delta==1)/.bound(((1-validmat$gHat1W)*validmat$dbar0W),length(which(data$S %in% comparisons[[s]])), bounds)
      H.AW <- as.numeric(data$A==1 & data$Delta==1) - as.numeric(data$A==0 & data$Delta==1)
      H.1W <- rep(1, nrow(data))
      H.0W <- rep(-1, nrow(data))

    } else{
      wt <- rep(1, nrow(data))
      H.AW <- as.numeric(data$A==1 & data$Delta==1)/.bound((validmat$gHat1W*validmat$dbar1W),length(which(data$S %in% comparisons[[s]])), bounds) - as.numeric(data$A==0 & data$Delta==1)/.bound(((1-validmat$gHat1W)*validmat$dbar0W),length(which(data$S %in% comparisons[[s]])), bounds)
      H.1W <- 1/.bound((validmat$gHat1W*validmat$dbar1W),length(which(data$S %in% comparisons[[s]])), bounds)
      H.0W <- -1/.bound(((1-validmat$gHat1W)*validmat$dbar0W),length(which(data$S %in% comparisons[[s]])), bounds)

    }


    if(fluctuation == "logistic"){
      logitUpdate<- glm(validmat$Yscale[which(data$Delta==1 & (data$S %in% comparisons[[s]]))] ~ -1 + offset(qlogis(validmat$QbarAW[which(data$Delta==1 & (data$S %in% comparisons[[s]]))])) +  H.AW[which(data$Delta==1 & (data$S %in% comparisons[[s]]))], family='quasibinomial', weights = wt[which(data$Delta==1 & (data$S %in% comparisons[[s]]))])

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

    out$clevercov[[s]] <- rep(NA, nrow(data))

    for(v in 1:V){
      out$psi[[v]][s] <- mean((validmat$Qbar1W.star - validmat$Qbar0W.star)[which(validmat$v==v & (data$S %in% comparisons[[s]]))])
      out$EICay[which(data$v==v & (data$S %in% comparisons[[s]])),(length(comparisons)*(v-1)+s)] <- ((wt*H.AW)*(data$Y - validmat$QbarAW.star) + validmat$Qbar1W.star - validmat$Qbar0W.star - out$psi[[v]][s])[which(data$v==v & (data$S %in% comparisons[[s]]))]/((length(which(data$v==v & (data$S %in% comparisons[[s]]))))/nrow(data))
      out$clevercov[[s]][which(data$v==v & (data$S %in% comparisons[[s]]) & data$Delta==1)] <- (wt*H.AW)[which(data$v==v & (data$S %in% comparisons[[s]]) & data$Delta==1)]
    }

    out$clevercov[[s]] <- stats::na.omit(out$clevercov[[s]])

    if(s==1){
      pooledVar <- list()
      for(v in 1:V){
        pooledVar[[v]] <- var(((wt*H.AW)*(data$Y - validmat$QbarAW.star) + validmat$Qbar1W.star - validmat$Qbar0W.star - out$psi[[v]][s])[which(data$v==v & (data$S %in% comparisons[[s]]))])/((length(which(data$S %in% comparisons[[s]]))))
      }
      out$Var <- mean(unlist(pooledVar))
    }


  }

  return(out)
}

#function for sampling from the estimated limit distribution

limitdist_sample <- function(V, bvt, NCO, EICpsipound, EICnco, var_ay, limitdist, data, comparisons){
  out <- list()
  psipoundvec <- NA
  for(v in 1:V){
    psipoundvec <- c(psipoundvec,bvt[[v]]$bias)
  }
  psipoundvec <- psipoundvec[-1]

  if(is.null(NCO)==FALSE){
    psipoundplusphivec <- NA
    for(v in 1:V){
      psipoundplusphivec <- c(psipoundplusphivec,(bvt[[v]]$bias + bvt[[v]]$bias_nco))
    }
    psipoundplusphivec <- psipoundplusphivec[-1]

    #overall covariance matrix for ztilde_poundplusphi
    EICpoundplusphi <- EICpsipound+EICnco
    EICmat_poundplusphi <- cbind(EICpoundplusphi, limitdist$EICay)
    out$covMat_poundplusphi <- (t(EICmat_poundplusphi)%*%EICmat_poundplusphi)/nrow(data)

    ztilde_poundplusphi_samp <- mvrnorm(n = 1000, mu=rep(0,ncol(EICmat_poundplusphi)), Sigma=out$covMat_poundplusphi/nrow(data))
  }

  #overall covariance matrix for ztilde
  EICmat <- cbind(EICpsipound, limitdist$EICay)

  out$covMat <- (t(EICmat)%*%EICmat)/nrow(data)

  #sample from multivariate ztilde
  ztilde_samp <- mvrnorm(n = 1000, mu=rep(0,ncol(EICmat)), Sigma=out$covMat/nrow(data))

  #selector for each sample
  biassample_psipound <- ztilde_samp[,(1:as.numeric(length(comparisons)*V))]
  if(is.null(NCO)==FALSE){
    biassample_psipoundplusphi <- ztilde_poundplusphi_samp[,(1:as.numeric(length(comparisons)*V))]
  }

  lambdatildeb2v <- matrix(NA, nrow=1000, ncol=length(psipoundvec))
  if(is.null(NCO)==FALSE){
    lambdatildencobias <- matrix(NA, nrow=1000, ncol=length(psipoundplusphivec))
  }
  for(b in 1:1000){
    lambdatildeb2v[b,] <- (biassample_psipound[b,]+psipoundvec)^2 + var_ay
    if(is.null(NCO)==FALSE){
      lambdatildencobias[b,] <- (biassample_psipoundplusphi[b,] + psipoundplusphivec)^2 + var_ay
    }
  }

  psisamp <- ztilde_samp[,(((as.numeric(length(comparisons)*V)+1)):(2*as.numeric(length(comparisons)*V)))]
  if(is.null(NCO)==FALSE){
    psisamp_poundplusphi <- ztilde_poundplusphi_samp[,(((as.numeric(length(comparisons)*V)+1)):(2*as.numeric(length(comparisons)*V)))]
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
  out$psi_pstarnv_b2v <- vector()
  psi_pstarnv_b2v_v <- list()
  out$psi_pstarnv_nco <- vector()
  psi_pstarnv_nco_v <- list()
  for(b in 1:1000){
    psi_pstarnv_b2v_v[[b]] <- vector()
    psi_pstarnv_nco_v[[b]] <- vector()
    for(v in 1:V){
      psi_pstarnv_b2v_v[[b]][v] <- sample_psi_pstarnv[[b]][v,which(lambdatildeb2v[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]==min(lambdatildeb2v[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]))]
      if(is.null(NCO)==FALSE){
        psi_pstarnv_nco_v[[b]][v] <- sample_psi_pstarnv[[b]][v,which(lambdatildencobias[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]==min(lambdatildencobias[b,((length(comparisons)*(v-1)+1):(length(comparisons)*(v)))]))]
      }
    }
    out$psi_pstarnv_b2v[b] <- mean(psi_pstarnv_b2v_v[[b]])
    if(is.null(NCO)==FALSE){
      out$psi_pstarnv_nco[b] <- mean(psi_pstarnv_nco_v[[b]])
    }
  }
  return(out)
}

