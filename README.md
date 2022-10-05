# EScvtmle

This package implements the experiment-selector cross-validated targeted maximum likelihood estimator (CV-TMLE) for integrating observational and RCT data described in Dang et al. (2022) [CITE]. The ES-CVTMLE selects and then analyzes the experiment (RCT only or RCT with real-world data) with the best estimated bias-variance tradeoff. If a negative control outcome (NCO) is available, the bias estimate may include the estimated ATE of treatment on the NCO, which facilitates inclusion of unbiased real-world data.
