README: Test Datasets for ACIC 2019

Files for the Data Challenge will be in the same form as the sample files here that have filenames  "testdatasetX.csv". 

The data are in the form (Y, A, V1, ... Vp), where
Y is the outcome (either binary or continuous) 
A is a binary treatment indicator
V1 through Vp are covariates. (The number of covariates, p, varies across datasets)


Files with the name "testdatasetX_cf.csv" have three columns. 
ATE = the population average treatment effect (PATE), defined as E(EY(1) - EY(0)).  The expectation is with respect to the population distribution of covariates. The population level effect is not identical to the sample average treatment effect.
EY1_i = expected value of the counterfactual outcome under treatment, conditional on covariates
EY0_i = expected value of the counterfactual outcome under no treatment, conditional on covariates

Note:  
When the observed Y is continuous, Y = A*EY1_i + (1-A) * EY0_i + epsilon, where epsilon is an error term having mean 0.  

When the observed Y is binary, Y is generated as a Bernoulli random variable, with probability  = A*EY1_i + (1-A) * EY0_i 