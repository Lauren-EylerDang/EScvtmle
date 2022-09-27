#'
#bound denominator of clever covariates
.bound <- function(x,n){
  if(any(x==0)){
    stop("Denominator of clever covariate 0 for at least one observation. Positivity violation?")
  }
  x <- pmax((5/(n)^(1/2)/log(n)), pmin(1,x))
  return(x)
}

#' @importFrom dplyr rename
# Function to pre-process data
# Includes removal of observations from observational dataset with W covariates not represented in RCT if txinrwd==FALSE
preprocess <- function(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, adjustnco = adjustnco){

  #remove observations missing treatment
  data <- rename(data, A = treatment_var)

  if(length(which(is.na(data$A)))>0) message("Removing observations with missing treatment variable.")
  data <- data[which(is.na(data$A)==FALSE),]

  #make A coded as 1=treatment of interest, 0=control
  data$A <- ifelse(data$A == treatment, 1, 0)

  data <- rename(data, S = study)
  data <- rename(data, Y = outcome)

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
    data <- rename(data, nco = NCO)

    if(adjustnco == TRUE & txinrwd == FALSE){
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

  return(data)
}


#apply selector_func to different datasets
apply_selector_func <- function(txinrwd, train, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO=NULL, Delta=NULL, Delta_NCO=NULL, adjustnco=adjustnco, target.gwt=target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL, comparisons){
  out <- list()
  for(s in 1:(length(comparisons))){

    #train regressions on A=0 only
    train_s <- train[which(train$S %in% comparisons[[s]]),]
    train_s$S[which(train_s$S!=1)]<-0

    if(txinrwd==TRUE){
      out[[s]] <- selector_func_txrwd(train_s = train_s, data=data, Q.SL.library, d.SL.library, g.SL.library, pRCT = pRCT, family = family, family_nco = family_nco, fluctuation = fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL)
    } else {
      out[[s]] <- selector_func_notxrwd(train_s = train_s, data=data, Q.SL.library, d.SL.library, g.SL.library, pRCT = pRCT, family = family, family_nco = family_nco, fluctuation = fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL=Q.discreteSL, d.discreteSL=d.discreteSL, g.discreteSL=g.discreteSL)
    }
  }
  return(out)
}

#' @importFrom stats predict
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

