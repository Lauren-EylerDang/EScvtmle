#' @title ES.cvtmle
#'
#' @description This function runs the experiment-selector cross-validated targeted maximum likelihood estimator (ES-CVTMLE) (Dang et al. 2022) for selecting and analyzing an optimal experiment, where candidate experiments include a randomized controlled trial (RCT) with or without various real-world datasets (RWD).
#'
#' @param txinrwd Whether active treatment is available in RWD (TRUE/FALSE). If FALSE, only the control arm of the RCT will be augmented with external data.
#' @param data The dataset
#' @param study Character name of variable indicating study participation (e.g. "S"). This variable should take a value of 1 for the RCT and should take a value of 0 for the RWD. Note that the code is currently set up only to handle two studies, but may be expanded to handle multiple studies in the future.
#' @param covariates Vector of character names of covariates to be adjusted for (e.g. c("W1", "W2"))
#' @param treatment_var Character name of treatment variable (e.g. "A")
#' @param treatment Value of treatment variable that corresponds to the active treatment (e.g. "DrugName" or 1). All other values of the treatment variable are assumed to be control.
#' @param outcome Character name of outcome variable (e.g. "Y"). If the outcome is a binary variable subject to censoring, censored observations should either be coded as NA or should be coded as 0 and a missingness indicator should be included (see parameter Delta below).
#' @param NCO Character name of negative control outcome variable (e.g. "nco") or NULL if no NCO available. If the NCO is a binary variable subject to censoring, censored observations should either be coded as NA or should be coded as 0 and a missingness indicator should be included (see parameter Delta_NCO below).
#' @param Delta Character name of a variable that is 0 if an observation was censored (missing outcome) and 1 otherwise. Missing outcomes may also be coded as NA, in which case a Delta variable will be added internally. If no missing outcomes, set Delta=NULL.
#' @param Delta_NCO Character name of a variable that is 0 if the value of NCO is missing and 1 otherwise. Missing NCOs may also be coded as NA, in which case a Delta_NCO variable will be added internally. If no missing NCO or no NCO, set Delta_NCO=NULL.
#' @param id ID variable for the independent unit
#' @param pRCT The probability of randomization to treatment in the RCT
#' @param V Number of cross-validation folds (default 10)
#' @param Q.SL.library Candidate algorithms for SuperLearner estimation of outcome regressions
#' @param d.SL.library.RCT Candidate algorithms for SuperLearner estimation of missingness mechanism for RCT-only
#' @param d.SL.library.RWD Candidate algorithms for SuperLearner estimation of missingness mechanism for RCT+RWD
#' @param g.SL.library Candidate algorithms for SuperLearner estimation of treatment mechanism for combined RCT/RWD analysis
#' @param Q.discreteSL Should a discrete SuperLearner be used for estimation of outcome regressions? (TRUE/FALSE)
#' @param d.discreteSL Should a discrete SuperLearner be used for estimation of missingness mechanism? (TRUE/FALSE)
#' @param g.discreteSL Should a discrete SuperLearner be used for estimation of treatment mechanism? (TRUE/FALSE)
#' @param family Either "binomial" for binary outcomes or "gaussian" for continuous outcomes
#' @param family_nco Family for negative control outcome
#' @param fluctuation 'logistic' (default for binary and continuous outcomes), or 'linear' describing fluctuation for targeted maximum likelihood estimation (TMLE) updating. If 'logistic' with a continuous outcome, outcomes are scaled to (0,1) for TMLE targeting and then returned to the original scale for parameter estimation.
#' @param comparisons A list of the values of the study variable that you would like to compare. For example, if you have an RCT labeled S=1 and RWD labeled S=0, you would use comparisons = list(c(1),c(1,0)) to compare RCT only to RCT + RWD. The first element of comparisons must be c(1) for the RCT only.
#' @param adjustnco Should we adjust for the NCO as a proxy of bias in the estimation of the ATE of A on Y? (TRUE/FALSE). Default is FALSE.
#' @param target.gwt As in the tmle R package (Gruber & van der Laan, 2012), if target.gwt is TRUE, the treatment mechanism is moved from the denominator of the clever covariate to the weight when fitting the coefficient for TMLE updating. Default TRUE.
#' @param bounds Optional bounds for truncation of the denominator of the clever covariate. The default is c(5/sqrt(n)/log(n),1).
#' @param cvControl A list of parameters to control the cross-validation process for the SuperLearners. See ?SuperLearner for more details.
#' @param MCsamp Number of Monte Carlo samples from the estimated limit distribution to use to estimate quantile-based confidence intervals. Default 1000.
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
#' @return Returns an object of class "EScvtmle" that is a list with the following components.
#' \describe{
#'  \item{ATE}{Average treatment effect (ATE) point estimates for the ES-CVTMLE estimator using the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#'  \item{foldATEs}{Average treatment effect (ATE) point estimates for each cross-validation fold of the ES-CVTMLE estimator using the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#'  \item{g}{g is a list of the same length as comparisons where each element of the list is a vector of the denominator of the covariate in front of the residual in the efficient influence curve for all observations in the experiment described by that element of comparisons. Values of g close to 0 or 1 indicate practical near-positivity violations.}
#'  \item{CI}{Estimated 95\% confidence intervals for the average treatment effect estimates of the ES-CVTMLE estimator using the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#'  \item{limitdistributionsample}{Monte Carlo samples for the average treatment effect estimates of the ES-CVTMLE estimator that are used to construct confidence intervals for the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#'  \item{Var}{Estimated variance of the ES-CVTMLE average treatment effect estimator using the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#'  \item{selected_byfold}{Vector noting which experiment from the list of comparisons was selected in each cross-validation fold of the ES-CVTMLE estimator using the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#'  \item{proportionselected}{Proportion of all cross-validation folds in which real-world (external) data were included in the analysis for the ES-CVTMLE estimator using the estimated bias squared plus variance selector ("b2v") and for the selector that includes an estimate of the ATE on a negative control outcome (NCO) in the bias term of the selector ("ncobias") if an NCO is available.}
#' }
#'
#' @details The experiment selector cross-validated targeted maximum likelihood estimator (ES-CVTMLE) aims to select the experiment that optimizes the bias-variance tradeoff for estimating a causal average treatment effect where different experiments may include a randomized controlled trial (RCT) alone or an RCT combined with real-world data.
#' Using cross-validation, the ES-CVTMLE separates the selection of the optimal experiment from the estimation of the ATE for the chosen experiment.
#' In order to avoid positivity violations, the package internally trims RWD so that no baseline covariate values are not represented in the RCT if active treatment is not available in the RWD.
#' The estimated bias term in the selector is a function of the difference in conditional mean outcome under control for the RCT compared to the combined experiment.
#' In order to help include truly unbiased external data in the analysis, the estimated average treatment effect on a negative control outcome may be added to the bias term in the selector by setting the parameter NCO to the character name of a negative control variable in the dataset.
#' For more details about this method, please see Dang et al. (2022).
#'
#' References:
#'
#' Dang LE, Tarp JM, Abrahamsen TJ, Kvist K, Buse JB, Petersen M, van der Laan M (2022). A Cross-Validated Targeted Maximum Likelihood Estimator for Data-Adaptive Experiment Selection Applied to the Augmentation of RCT Control Arms with External Data. arXiv:2210.05802 [stat.ME]
#'
#' Susan Gruber, Mark J. van der Laan (2012). tmle: An R Package for Targeted Maximum Likelihood Estimation. Journal of Statistical Software, 51(13), 1-35. URL <http://www.jstatsoft.org/v51/i13/>.
#'
#' @examples
#' \donttest{data(wash)
#' #For unbiased external controls, use:
#' dat <- wash[which(wash$study %in% c(1,2)),]
#' dat$study[which(dat$study==2)]<-0
#' set.seed(2022)
#' results_rwd1 <- ES.cvtmle(txinrwd=TRUE,
#'                           data=dat, study="study",
#'                           covariates=c("aged", "sex", "momedu", "hfiacat"),
#'                           treatment_var="intervention", treatment=1,
#'                           outcome="laz", NCO="Nlt18scale",
#'                           Delta=NULL, Delta_NCO=NULL,
#'                           pRCT=0.5, V=5, Q.SL.library=c("SL.glm"),
#'                           g.SL.library=c("SL.glm"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
#'                           family="gaussian", family_nco="gaussian", fluctuation = "logistic",
#'                           comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)
#' print.EScvtmle(results_rwd1)
#' }
#'
#' @export

ES.cvtmle <- function(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO=NULL, Delta=NULL, Delta_NCO=NULL, id=NULL, pRCT, V=10, Q.SL.library, d.SL.library.RCT, d.SL.library.RWD, g.SL.library, Q.discreteSL, d.discreteSL, g.discreteSL, family, family_nco, fluctuation = "logistic", comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE, bounds=NULL, cvControl = list(), MCsamp=1000){

  if (length(comparisons)>2) stop("Package currently compares two experiments. Check back for updates to compare multiple experiments.")

  if (any((data$study %in% c(0,1))==FALSE)) stop("Package currently considers two studies, where study=1 is an RCT, and study=0 is a real-world dataset.")

  if (comparisons[[1]]!=1) stop("First comparison should be c(1) (ie compare to RCT only).")

  data <- preprocess(txinrwd, data, study, covariates, treatment_var, treatment, outcome, NCO, Delta, Delta_NCO, id, adjustnco)
  if(any(data$Delta==0)){
    Delta = "Delta"
  }

  if(any(data$NCO_delta==0)){
    Delta_NCO = "NCO_delta"
  }

  #Create cross-validation folds that preserve proportion of RCT (if txinwrd=TRUE) or of RCT controls (if txinrwd=FALSE) as well as proportion of outcomes in validation sets

    ids <- data$S

    if(txinrwd==TRUE){
      if(family=="binomial"){
        ids[which(data$S==1 & data$Y==1)] <- 0
      } else {
        ids[which(data$S==1)] <- 0
      }
    } else {
      if(family=="binomial"){
        ids[which(data$S==1 & data$A==0 & data$Y==1)] <- 0
      } else {
        ids[which(data$S==1 & data$A==0)] <- 0
      }
    }

    folds <- make_folds(data, fold_fun = folds_vfold, V=V, cluster_ids = data$id, strata_ids = ids)

  data$v <- rep(NA, nrow(data))
  for(v in 1:V){
    data$v[folds[[v]]$validation_set]<-v
  }

  #number of independent units
  n.id <- length(unique(data$id))

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
  EICpsipound <- matrix(0, nrow=n.id, ncol=length(comparisons)*V)
  EICnco <- matrix(0, nrow=n.id, ncol=length(comparisons)*V)

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

    selector[[v]] <- apply_selector_func(txinrwd, train, data, Q.SL.library, d.SL.library.RCT, d.SL.library.RWD, g.SL.library, pRCT, family, family_nco, fluctuation, NCO, Delta, Delta_NCO, adjustnco, target.gwt, Q.discreteSL, d.discreteSL, g.discreteSL, comparisons, bounds, cvControl)

    if(txinrwd==TRUE){
      bvt[[v]] <- bvt_txinrwd(v, selector, NCO, comparisons, train, data, fluctuation, family, n.id)
    } else {
      bvt[[v]] <- bvt_notxinrwd(v, selector, NCO, comparisons, train, data, fluctuation, family, n.id)
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
  limitdist <- limitdistvar(V, valid_initial, data, folds, family, fluctuation, Delta, pRCT, target.gwt, comparisons, bounds, n.id)

  results$g <- list()
  for(s in 1:length(comparisons)){
    results$g[[s]] <- as.vector(1/(abs(limitdist$clevercov[[s]])))
  }


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
  } else {
    results$ATE$ncobias <- results$foldATEs$ncobias <- NULL
  }

  #sample from limit distribution
  limitdistsamp <- limitdist_sample(V, bvt, NCO, EICpsipound, EICnco, var_ay, limitdist, n.id, comparisons, MCsamp)

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
  } else {
    results$Var$ncobias <- results$CI$ncobias <- NULL
  }

  results$selected_byfold <- list()
  results$selected_byfold$b2v <- unlist(proportionselected$b2v)
  results$selected_byfold$ncobias <- unlist(proportionselected$ncobias)

  results$proportionselected <- list()
  results$proportionselected$b2v <- (mean(unlist(proportionselected$b2v))-1)
  if(is.null(NCO)==FALSE){
    results$proportionselected$ncobias <- (mean(unlist(proportionselected$ncobias))-1)
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
#' @details Prints the average treatment effect (ATE) point estimate and 95\% confidence interval for the ES-CVTMLE estimator (object of class "EScvtmle") using the estimated bias squared plus variance experiment-selector. If a negative control outcome (NCO) is available, this function also prints the ATE point estimate and 95\% confidence interval for the selector that includes the estimated ATE on the NCO in the bias term. See Dang et al. (2022) <arXiv:2210.05802> for more details.
#' @return No return value. Called to print a summary of the results for objects of class "EScvtmle".
#' @export print.EScvtmle
#' @export
print.EScvtmle <- function(x, ...) {
  if(identical(class(x), "EScvtmle")){
    if(is.null(x$NCO)==FALSE){
      cat("Experiment-Selector CV-TMLE Average Treatment Effect Estimate")
      cat("\n   Without NCO: ", paste(round(x$ATE$b2v, 3), " 95% CI (", round(x$CI$b2v[1],3), " - ", round(x$CI$b2v[2],3), ")", sep=""))
      cat("\n   RWD included in ", paste((x$proportionselected$b2v*100), "% of folds.", sep=""),"\n")

      cat("\n Experiment-Selector CV-TMLE Average Treatment Effect Estimate")
      cat("\n   With NCO: ", paste(round(x$ATE$ncobias, 3), " 95% CI (", round(x$CI$ncobias[1],3), " - ", round(x$CI$ncobias[2],3), ")", sep=""))
      cat("\n   RWD included in ", paste((x$proportionselected$ncobias*100), "% of folds.", sep=""),"\n")
    } else {
      cat("Experiment-Selector CV-TMLE Average Treatment Effect Estimate")
      cat("\n   Without NCO: ", paste(round(x$ATE$b2v, 3), " 95% CI (", round(x$CI$b2v[1],3), " - ", round(x$CI$b2v[2],3), ")", sep=""))
      cat("\n   RWD included in ", paste((x$proportionselected$b2v*100), "% of folds.", sep=""),"\n")
    }
  } else {
    stop("Error: Object class is not EScvtmle \n")
  }
}

#' @title plot.EScvtmle
#'
#' @description Plots fold-specific average treatment effect (ATE) estimates and a histogram of Monte Carlo sample ATE estimates used to construct confidence intervals.
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
#' @return Returns a graphical object of class "grob" that contains two side-by-side plots: one of the fold-specific average treatment effect estimates for all cross-validation folds (including information regarding which experiment was selected in each fold), and the other of a histogram of the Monte Carlo samples that are used to construct confidence intervals. If a negative control outcome (NCO) is available, the plots are for the selector that includes the estimated average treatment effect on the NCO in the bias estimate. If not, the plots are for the selector that uses the estimated bias squared plus the variance selector, without information from an NCO. For more information about the different selectors, the use of cross-validation, or the construction of confidence intervals for this method, please see Dang et al. (2022)  <arXiv:2210.05802>.
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
        "Experiment"=as.factor(x$selected_byfold$b2v)
      )
      xdf <- data.frame("Samples"=x$limitdistributionsample$b2v)
      plot1 <- ggplot(df, aes(x=Fold, y=ATE, shape=Experiment, color=Experiment)) + geom_point() + geom_hline(yintercept=x$ATE$b2v) + ggtitle("ATE Estimates by Fold") + theme(plot.title = element_text(hjust = 0.5))
      plot2 <- ggplot(xdf, aes(x=Samples)) + geom_histogram() + ggtitle("Histogram of Monte Carlo Samples for Confidence Interval Estimation") +labs(y= "Frequency", x = "ATE Estimate", caption = "Red Lines Mark 95% CI") + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = x$CI$b2v, color="red")

      grid.arrange(plot1, plot2, ncol=2)
    }
  } else {
    stop("Error: Object class is not EScvtmle \n")
  }
}
