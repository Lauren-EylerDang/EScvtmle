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
  A[which(S==2)] <- rbinom(n2, 1, pRCT)

  #Outcome regression

  Y <- rbinom(N, 1, plogis(0.1*W1 - 0.5*W2 +bA*A))
  nco <- rbinom(N, 1, plogis(0.5*W1 - 0.1*W2))

  data <- data.frame(S, W1, W2, A, Y, nco)

  return(data)
}

data <- generateData(n1, n2, 0.5)

library(SuperLearner)

results <- ES.cvtmle(txinrwd=TRUE,data=data, study="S",
                     covariates=c("W1", "W2"),
                     treatment_var="A", treatment=1,
                     outcome="Y", NCO="nco",
                     Delta=NULL, Delta_NCO=NULL,
                     pRCT=0.5, V=10, Q.SL.library=c("SL.glm"),
                     g.SL.library=c("SL.glm"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
                     family="binomial", family_nco="binomial", fluctuation = "logistic",
                     comparisons = list(c(1),c(1,2)), adjustnco = FALSE, target.gwt = TRUE)

test_that("Correct number of folds", {
  expect_equal(length(results$foldATEs$b2v), 10)
})

test_that("Correct selection of ATE by experiment across folds", {
  expect_equal(results$foldATEs$b2v[which(results$selected_byfold$b2v==results$selected_byfold$ncobias)], results$foldATEs$ncobias[which(results$selected_byfold$b2v==results$selected_byfold$ncobias)])
})

test_that("Output for which experiment selected in each fold correct", {
  expect_true(all(results$selected_byfold$b2v %in% c(1,2)))
  expect_true(all(results$selected_byfold$ncobias %in% c(1,2)))
})


