# On the estimation of population size – A comparison of capture-recapture and multiplier-benchmark methods

This directory contains the code to produce the simulation work for the thesis project "On the estimation of population size – A comparison of capture-recapture and multiplier-benchmark methods". It contains code to run the simulations comparing the capture-recapture approach to the multiplier-benchmark approach in multiple scenarios, including both the assumptions met and violated situations. 

## Framework

### FrameworkCRCvsMBM-Simu1.Rmd

This markdown file contains the codes comparing the sampling mechanisms between the capture-recapture and the multiplier-benchmark methods, including the simulation goals, setup, simulation codes, and codes for the plots.

### Simu1.Scenarios.Functions.R

This program contains functions computing the estimates of the total size $N$, the relative bias and variances based on different estimators in the two methods, as well as the visualization codes for the comparison.

## Misuse 1

### Misuse-Simu2.Rmd

This markdown file contains the codes exploring the multiplier-benchmark method regarding the impact of the samples coverage on the accuracy and precision of the estimated $N$ when the benchmark sample has undercounting issue. The investigation is done under different level of the sources dependence. 

### Simu2.Scenario.Function.R

This program contains corresponding functions computing the estimates of the total size $N$, the relative bias and variances based on different estimators in the two methods, as well as the visualization codes for the comparison.

## Misuse 2

### Misuse-Simu3.Rmd

This markdown file contains the codes exploring the capture-recapture method regarding the impact of the sample coverage on the accuracy and precision of the estimated $N$ when some subjects are structurally missed in the sampling process for one of the samples.

### Simu3.Scenario.Function.R

This program contains corresponding functions computing the estimates of the total size $N$, the relative bias and variances based on different estimators in the two methods, as well as the visualization codes for the comparison.

