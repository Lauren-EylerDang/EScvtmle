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
#' @importFrom dplyr rename
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
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
#' print.EScvtmle(results_rwd1)
#' @export

ES.cvtmle <- function(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, pRCT, V=10, Q.SL.library, d.SL.library, g.SL.library, Q.discreteSL, d.discreteSL, g.discreteSL, family, family_nco, fluctuation = "logistic", comparisons = list(c(1),c(1,2)), adjustnco = FALSE, target.gwt = TRUE){

  if (length(comparisons)>2) stop("Package currently compares two experiments. Check back for updates to compare multiple experiments.")

  if (comparisons[[1]]!=1) stop("First comparison should be c(1) (ie compare to RCT only).")

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
      var_ay[(length(comparisons)*(v-1)+s)] <- bvt[[v]]$var[s]
    }


  }

  valid_initial <- validpreds(data, folds, V, selector, pRCT, Delta, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons)

  results$ATE <- list()
  results$ATE$b2v <- vector()
  results$ATE$ncobias <- vector()

  results$foldATEs <- list()
  results$foldATEs$b2v <- vector()
  results$foldATEs$ncobias <- vector()

  #estimate components of limit distribution
  limitdist <- limitdistvar(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons)

  pool <- vector()
  for(v in 1:V){
    pool[v]<- limitdist$psi[[v]][proportionselected$b2v[[v]]]
  }
  results$foldATEs$b2v <- pool
  results$ATE$b2v <- mean(pool)

  if(is.null(NCO)==FALSE){
    pool <- vector()
    for(v in 1:V){
      pool[v]<- limitdist$psi[[v]][proportionselected$ncobias[[v]]]
    }
    results$foldATEs$ncobias <- pool
    results$ATE$ncobias <- mean(pool)
  }

  #sample from limit distribution
  limitdistsamp <- limitdist_sample(V, bvt, NCO, EICpsipound, EICnco, var_ay, limitdist, data, comparisons)

  results$CI$b2v <- list()
  results$CI$ncobias <- list()

  results$limitdistributionsample <- list()
  results$limitdistributionsample$b2v <- results$ATE$b2v + limitdistsamp$psi_pstarnv_b2v
  results$limitdistributionsample$nco <- results$ATE$ncobias + limitdistsamp$psi_pstarnv_nco

  #use quantiles of samples from limit distribution to estimate CI unless only RCT selected in all folds, then use standard CV-TMLE IC-based variance estimates
  if(any(unlist(proportionselected$b2v)!=1)){
    results$Var$b2v <- var(limitdistsamp$psi_pstarnv_b2v)
    results$CI$b2v <- results$ATE$b2v + quantile(limitdistsamp$psi_pstarnv_b2v, probs = c(0.025,0.975))
  } else {
    results$Var$b2v <- limitdist$Var
    results$CI$b2v <- c((results$ATE$b2v - 1.96*(limitdist$Var)^(1/2)), (results$ATE$b2v + 1.96*(limitdist$Var)^(1/2)))
  }

  if(is.null(NCO)==FALSE){
    if(any(unlist(proportionselected$ncobias)!=1)){
      results$Var$ncobias <- var(limitdistsamp$psi_pstarnv_nco)
      results$CI$ncobias <- results$ATE$ncobias + quantile(limitdistsamp$psi_pstarnv_nco, probs = c(0.025,0.975))
    } else {
      results$Var$ncobias <- limitdist$Var
      results$CI$ncobias <- c((results$ATE$ncobias - 1.96*(limitdist$Var)^(1/2)), (results$ATE$ncobias + 1.96*(limitdist$Var)^(1/2)))
    }
  }

  results$selected_byfold <- list()
  results$selected_byfold$b2v <- unlist(proportionselected$b2v)
  results$selected_byfold$ncobias <- unlist(proportionselected$ncobias)

  results$proportionselected_mean <- list()
  results$proportionselected_mean$b2v <- (mean(unlist(proportionselected$b2v))-1)
  if(is.null(NCO)==FALSE){
    results$proportionselected_mean$ncobias <- (mean(unlist(proportionselected$ncobias))-1)
  }
  results$NCO <- NCO

  class(results) <- "EScvtmle"

  return(results)
}

#' @title print.EScvtmle
#'
#' @description Prints output from object produced by ES.cvtmle function
#'
#' @param x An object of class "EScvtmle"
#' @param ... Other arguments to print
#' @method print EScvtmle
#' @export print.EScvtmle
#' @export
print.EScvtmle <- function(x, ...) {
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

#' @title plot.EScvtmle
#'
#' @description Plots fold-specific ATE estimates and histogram of monte carlo sample ATE estimates
#'
#' @param x An object of class "EScvtmle"
#' @param ... Other arguments to plot
#' @method plot EScvtmle
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom gridExtra grid.arrange
#' @export plot.EScvtmle
#' @export
plot.EScvtmle <- function(x, ...) {
  if(identical(class(x), "EScvtmle")){
    Fold <- ATE <- Experiment <- Samples <- NULL
    if(is.null(x$NCO)==FALSE){
      df <- data.frame(
        "Fold"=seq(1,length(x$selected_byfold$ncobias),1),
        "ATE"=x$foldATEs$ncobias,
        "Experiment"=as.factor(x$selected_byfold$ncobias)
      )
      xdf <- data.frame("Samples"=x$limitdistributionsample$nco)
      plot1 <- ggplot(df, aes(x=Fold, y=ATE, shape=Experiment, color=Experiment)) + geom_point() + geom_hline(yintercept=x$ATE$ncobias) + ggtitle("ATE Estimates by Fold") + theme(plot.title = element_text(hjust = 0.5))
      plot2 <- ggplot(xdf, aes(x=Samples)) + geom_histogram() + ggtitle("Histogram of Monte Carlo Samples") +labs(y= "Frequency", x = "ATE Estimate", caption = "Red Lines Mark 95% CI") + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = x$CI$ncobias, color="red")

      grid.arrange(plot1, plot2, ncol=2)

    } else {
      df <- data.frame(
        "Fold"=seq(1,length(x$selected_byfold$b2v),1),
        "ATE"=x$foldATEs$b2v,
        "Experiment"=x$selected_byfold$b2v
      )
      xdf <- data.frame("Samples"=x$limitdistributionsample$b2v)
      plot1 <- ggplot(df, aes(x=Fold, y=ATE, shape=Experiment, color=Experiment)) + geom_point() + geom_hline(yintercept=x$ATE$b2v) + ggtitle("ATE Estimates by Fold") + theme(plot.title = element_text(hjust = 0.5))
      plot2 <- ggplot(xdf, aes(x=Samples)) + geom_histogram() + ggtitle("Histogram of Monte Carlo Samples") +labs(y= "Frequency", x = "ATE Estimate", caption = "Red Lines Mark 95% CI") + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = x$CI$b2v, color="red")

    }
  } else {
    stop("Error: Object class is not EScvtmle \n")
  }
}
