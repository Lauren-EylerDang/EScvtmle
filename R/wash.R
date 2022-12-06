#' WASH Benefits Bangladesh Dataset
#'
#' This dataset was constructed from the publicly-available WASH Benefits Bangladesh cluster randomized controlled trial (RCT) dataset.
#' The results of this trial were originally reported by Luby et al. (2018), and the original dataset may be downloaded from https:://osf.io.wvyn4/.
#' The trial found no evidence of an effect of an intervention to improve sanitation, including construction of improved latrines, on child length-for-age Z-scores (laz).
#' A subsequent re-analysis by Arnold et al. (2018) of the control arm of this dataset as an observational cohort did find an effect of having an improved latrine at baseline on child laz.
#' The authors concluded that observational analyses of water, sanitation, and hygiene (WASH) interventions may suffer from unmeasured confounding.
#' To demonstrate how conducting a small RCT combined with unbiased or biased observational data could prevent unmeasured confounding from influencing results, we construct the dataset for this software package, as follows.
#' Study 1: A random sample of 150 "Sanitation" arm observations and 150 "Control" arm observations with complete information was taken from the overall RCT, with "study" variable set to 1.
#' Study 2: A second random sample of 150 "Sanitation" arm observations and 150 "Control" arm observations with complete information was taken from the remaining RCT observations, with "study" variable set to 2 to mimic an unbiased external dataset.
#' Study 3: From the "Control" arm observations not included in the study=1 sample, 150 observations who had improved latrines at baseline and 150 observations who did not have improved latrines at baseline were sampled, with "study" variable set to 3 to mimic a biased external dataset.
#' The data contained in this file consist of all three "studies".
#' Because this study was not set up to have a negative control outcome for length-for-age Z-score, the options were limited.
#' We would like a variable that is associated with socioeconomic status (SES) because that is a likely cause of the unmeasured confounding highlighted by Arnold et al. (2018).
#' We chose number of household members <=18 years old as an NCO, because prior studies have shown this variable to be associated with SES in Bangladesh (The World Bank, 2013),
#' but it is unlikely to be affected by having an improved latrine.
#' We scaled this variable to match the scale of the true outcome (length-for-age Z-score).
#'
#' @docType data
#'
#' @usage data(wash)
#'
#' @format An object of class "data.frame"
#' \describe{
#'  \item{intervention}{For studies 1 and 2: 1 if participant was in the "Sanitation" arm and 0 if participant was in the "Control" arm. For study 3: 1 if participant's household had an improved latrine at baseline and 0 otherwise.}
#'  \item{study}{Study variable indicating RCT sample or external dataset as described above.}
#'  \item{laz}{Child length-for-age Z-score at 2 years post-baseline.}
#'  \item{aged}{Child's age in days.}
#'  \item{sex}{Child's sex.}
#'  \item{momedu}{Mother's education level.}
#'  \item{hfiacat}{Category of household food insecurity. Levels are "Food Secure", "Mildly Food Insecure", "Moderately to Severely Food Insecure".}
#'  \item{Nlt18scale}{Scaled number of household members <= 18 years old.}
#' }
#' @references
#' Luby SP, Rahman M, Arnold BF, et al. Effects of water quality, sanitation, handwashing, and nutritional interventions on diarrhoea and child growth in rural Bangladesh: a cluster randomised controlled trial. The Lancet Global Health. 2018;6(3):e302-e315. doi:10.1016/S2214-109X(17)30490-4
#'
#' Arnold BF, Null C, Luby SP, Colford JM. Implications of WASH Benefits trials for water and sanitation – Authors’ reply. The Lancet Global Health. 2018;6(6):e616-e617. doi:10.1016/S2214-109X(18)30229-8
#'
#' The World Bank. (2013). Bangladesh Poverty Assessment: Assessing a Decade of Progress in Reducing Poverty, 2000-2010. Bangladesh Development Series. Paper No. 31. https://documents1.worldbank.org/curated/en/109051468203350011/pdf/785590NWP0Bang00Box0377348B0PUBLIC0.pdf
#'
#' @source \url{https://osf.io/wvyn4/}
#' @examples
#'
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
"wash"
