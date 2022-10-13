# EScvtmle

This package implements the experiment-selector cross-validated targeted maximum likelihood estimator (ES-CVTMLE) for integrating observational and RCT data described in Dang et al. (2022). The ES-CVTMLE selects and then analyzes the experiment (RCT only or RCT with real-world data) with the best estimated bias-variance tradeoff. If a negative control outcome (NCO) is available, the bias estimate may include the estimated ATE of treatment on the NCO, which facilitates inclusion of unbiased real-world data. For more information, see:

Dang LE, Tarp JM, Abrahamsen TJ, Kvist K, Buse JB, Petersen M, van der Laan M (2022). A Cross-Validated Targeted Maximum Likelihood Estimator for Data-Adaptive Experiment Selection Applied to the Augmentation of RCT Control Arms with External Data. arXiv:2210.05802 [stat.ME]
