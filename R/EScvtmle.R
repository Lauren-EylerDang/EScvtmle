#' @title ES.cvtmle
#'
#' @description Runs Experiment-Selector CV-TMLE for selecting and analyzing optimal experiment as RCT with or without RWD
#'
#' @param txinrwd Whether active treatment is available in RWD (TRUE/FALSE)
#' @param data The dataset
#' @param study Character name of variable indicating study participation (e.g. "S"). This variable should take values of 1 for the RCT and can take other values for the other study. Note that the code is currently set up only to handle two studies, but may be expanded to handle multiple studies in the future.
#' @param covariates Vector of character names of covariates to be adjusted for (e.g. c("W1", "W2"))
#' @param treatment_var Character name of treatment variable (e.g. "A")
#' @param treatment Name of treatment of interest (e.g. "DrugName" or 1)
#' @param outcome Character name of outcome variable (e.g. "Y"). If Y is a binary variable subject to censoring, it should be coded as 0 for observations that were censored.
#' @param NCO Character name of negative control outcome variable (e.g. "nco") or NULL if no NCO available. If NCO is a binary variable subject to censoring, it should be coded as 0 for observations that were censored.
#' @param Delta Character name of a variable that is 0 if an observation was censored (missing binary outcome) and 1 otherwise. Missing outcomes may also be coded as NA, in which case a Delta variable will be added internally. If no missing outcomes, set Delta=NULL.
#' @param Delta_NCO Character name of a variable that is 0 if the value of NCO is missing and 1 otherwise. Missing NCOs may also be coded as NA, in which case a Delta_NCO variable will be added internally. If no missing NCO or no NCO, set Delta_NCO=NULL.
#' @param pRCT The probability of randomization to treatment in the RCT
#' @param V Number of cross-validation folds (default 10).
#' @param Q.SL.library Candidate algorithms for SuperLearner estimation of outcome regressions
#' @param d.SL.library Candidate algorithms for SuperLearner estimation of missingness mechanism
#' @param g.SL.library Candidate algorithms for SuperLearner estimation of treatment mechanism for combined RCT/RWD analysis
#' @param Q.discreteSL Should a discrete SuperLearner be used for estimation of outcome regressions? (TRUE/FALSE)
#' @param d.discreteSL Should a discrete SuperLearner be used for estimation of missingness mechanism? (TRUE/FALSE)
#' @param g.discreteSL Should a discrete SuperLearner be used for estimation of treatment mechanism? (TRUE/FALSE)
#' @param family Either "binomial" for binary outcomes or "gaussian" for continuous outcomes
#' @param family_nco Family for negative control outcome
#' @param fluctuation 'logistic' (default for binary and continuos outcomes), or 'linear' describing fluctuation for TMLE updating. If 'logistic' with a continuous outcome, outcomes are scaled to (0,1) for TMLE targeting and then returned to the original scale for parameter estimation.
#' @param comparisons A vector of the values of study variable S that you would like to consider. For example, if you have an RCT labeled S=1 and RWD labeled S=2, you would use comparisons = list(c(1),c(1,2)) to compare RCT only to RCT + RWD.
#' @param adjustnco Should we adjust for the NCO as a proxy of bias in the estimation of the ATE of A on Y? (TRUE/FALSE)
#' @param target.gwt As in tmle package, when TRUE, move g from denominator of clever covariate to the weight when fitting coefficient for TMLE updating.
#'
#' @importFrom origami make_folds
#' @importFrom origami folds_vfold
#' @importFrom MASS mvrnorm
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom SuperLearner SuperLearner
#'
#' @return Returns the ATE estimate with 95% confidence intervals for the Experiment-Selector CV-TMLE and the proportion of folds in which RWD was included in the estimate.
#' @examples
#' data(wash)
#' #For unbiased external controls, use:
#' dat <- wash[which(wash$study %in% c(1,2)),]
#' library(SuperLearner)
#' set.seed(2022)
#' results_rwd1 <- ES.cvtmle(txinrwd=TRUE,
#'                           data=dat, study="study",
#'                           covariates=c("aged", "sex", "momedu", "hfiacat"),
#'                           treatment_var="treatment", treatment=1,
#'                           outcome="laz", NCO="Nlt18scale",
#'                           Delta=NULL, Delta_NCO=NULL,
#'                           pRCT=0.5, V=10, Q.SL.library=c("SL.glm"),
#'                           g.SL.library=c("SL.glm"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
#'                           family="gaussian", family_nco="gaussian", fluctuation = "logistic",
#'                           comparisons = list(c(1),c(1,2)), adjustnco = FALSE, target.gwt = TRUE)
#' print(results_rwd1)
#' @export

ES.cvtmle <- function(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, pRCT, V=10, Q.SL.library, d.SL.library, g.SL.library, Q.discreteSL, d.discreteSL, g.discreteSL, family, family_nco, fluctuation = "logistic", comparisons = list(c(1),c(1,2)), adjustnco = FALSE, target.gwt = TRUE){

  if (length(comparisons)>2) stop("Package currently compares two experiments. Check back for updates to compare multiple experiments.")

  if (comparisons[[1]]!=1) stop("First comparison should be c(1) (ie compare to RCT only).")

  data <- preprocess(data, study, covariates, treatment_var, treatment, outcome, NCO, Delta, Delta_NCO, adjustnco)

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

  EICay <- vector()
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

    selector[[v]] <- apply_selector_func(txinrwd, train, data, Q.SL.library, d.SL.library, g.SL.library, pRCT, family, family_nco, fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons)

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
      EICay[(length(comparisons)*(v-1)+s)] <- bvt[[v]]$var[s]
    }


  }

  valid_initial <- validpreds(data, folds, V, selector, pRCT, Delta, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons)

  results$ATE <- list()
  results$ATE$b2v <- vector()
  results$ATE$ncobias <- vector()

  limitdist <- limitdistvar(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons)

  pool <- vector()
  for(v in 1:V){
    pool[v]<- limitdist$psi[[v]][proportionselected$b2v[[v]]]
  }
  results$ATE$b2v <- mean(pool)

  if(is.null(NCO)==FALSE){
    pool <- vector()
    for(v in 1:V){
      pool[v]<- limitdist$psi[[v]][proportionselected$ncobias[[v]]]
    }
    results$ATE$ncobias <- mean(pool)
  }

  #Estimated covariance matrices
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
    covMat_poundplusphi <- (t(EICmat_poundplusphi)%*%EICmat_poundplusphi)/nrow(data)

    ztilde_poundplusphi_samp <- mvrnorm(n = 1000, mu=rep(0,ncol(EICmat_poundplusphi)), Sigma=covMat_poundplusphi/nrow(data))
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
  }

  lambdatildeb2v <- matrix(NA, nrow=1000, ncol=length(psipoundvec))
  if(is.null(NCO)==FALSE){
    lambdatildencobias <- matrix(NA, nrow=1000, ncol=length(psipoundplusphivec))
  }
  for(b in 1:1000){
    lambdatildeb2v[b,] <- (biassample_psipound[b,]+psipoundvec)^2 + EICay
    if(is.null(NCO)==FALSE){
      lambdatildencobias[b,] <- (biassample_psipoundplusphi[b,] + psipoundplusphivec)^2 + EICay
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
  psi_pstarnv_b2v <- vector()
  psi_pstarnv_b2v_v <- list()
  psi_pstarnv_nco <- vector()
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
    psi_pstarnv_b2v[b] <- mean(psi_pstarnv_b2v_v[[b]])
    if(is.null(NCO)==FALSE){
      psi_pstarnv_nco[b] <- mean(psi_pstarnv_nco_v[[b]])
    }
  }


  results$CI$b2v <- list()
  results$CI$ncobias <- list()

  if(any(unlist(proportionselected$b2v)!=1)){
    results$Var$b2v <- var(psi_pstarnv_b2v)
    results$CI$b2v <- results$ATE$b2v + quantile(psi_pstarnv_b2v, probs = c(0.025,0.975))
  } else {
    results$Var$b2v <- limitdist$Var
    results$CI$b2v <- c((results$ATE$b2v - 1.96*(limitdist$Var)^(1/2)), (results$ATE$b2v + 1.96*(limitdist$Var)^(1/2)))
  }

  if(is.null(NCO)==FALSE){
    if(any(unlist(proportionselected$ncobias)!=1)){
      results$Var$ncobias <- var(psi_pstarnv_nco)
      results$CI$ncobias <- results$ATE$ncobias + quantile(psi_pstarnv_nco, probs = c(0.025,0.975))
    } else {
      results$Var$ncobias <- limitdist$Var
      results$CI$ncobias <- c((results$ATE$ncobias - 1.96*(limitdist$Var)^(1/2)), (results$ATE$ncobias + 1.96*(limitdist$Var)^(1/2)))
    }
  }

  results$proportionselected_mean <- list()
  results$proportionselected_mean$b2v <- (mean(unlist(proportionselected$b2v))-1)
  if(is.null(NCO)==FALSE){
    results$proportionselected_mean$ncobias <- (mean(unlist(proportionselected$ncobias))-1)
  }
  results$NCO <- NCO

  class(results) <- "EScvtmle"

  return(results)
}

print.EScvtmle <- function(x,...) {
  if(identical(class(x), "EScvtmle")){
      if(is.null(x$NCO)==FALSE){
        cat("Experiment-Selector CV-TMLE Average Treatment Effect Estimate")
        cat("\n   Without NCO: ", paste(round(x$ATE$b2v, 3), " 95% CI (", round(x$CI$b2v[1],3), " - ", round(x$CI$b2v[2],3), ")", sep=""))
        cat("\n   RWD included in ", paste((x$proportionselected_mean$b2v*100), "% of folds.", sep=""),"\n")

        cat("\n Experiment-Selector CV-TMLE Average Treatment Effect Estimate")
        cat("\n   With NCO: ", paste(round(x$ATE$ncobias, 3), " 95% CI (", round(x$CI$ncobias[1],3), " - ", round(x$CI$ncobias[2],3), ")", sep=""))
        cat("\n   RWD included in ", paste((x$proportionselected_mean$ncobias*100), "% of folds.", sep=""),"\n")
      } else {
        cat("Experiment-Selector CV-TMLE Average Treatment Effect Estimate")
        cat("\n   Without NCO: ", paste(round(x$ATE$b2v, 3), " 95% CI (", round(x$CI$b2v[1],3), " - ", round(x$CI$b2v[2],3), ")", sep=""))
        cat("\n   RWD included in ", paste((x$proportionselected_mean$b2v*100), "% of folds.", sep=""),"\n")
      }
  } else {
    stop("Error: Object class is not EScvtmle \n")
  }
}
