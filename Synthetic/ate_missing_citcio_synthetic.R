# R script for simulation_ate_v10

args <- commandArgs(TRUE)
dir <- args[1]
N <- as.numeric(args[2]) # number of simulations
nb_cores <- as.numeric(args[3]) # set number of CPU cores 
missing <- args[4]
prob.miss <- as.numeric(args[5])
mim <- as.numeric(args[6])
Date <- args[7]

###################################################################################
# Load libraries
###################################################################################
library(MASS) # mvrnorm
library(caret) # cross-validated training
library(randomForest) # Random Forests
library(ranger) # Random Forests prediction
library(FactoMineR) # Factorial data analysis
library(grf) # Doubly robust treatment effect estimation
library(cobalt)

library(misaem) # Logistic regression with missing values
library(missMDA) # PCA/MCA with missing values + iterative PC imputation
library(missForest) # Random forests with missing values
library(mice) # Multiple imputation

library(foreach)
library(doParallel)


library(assertthat)
library(norm)
library(softImpute)

library(dplyr)
library(tidyr)
library(directlabels)
min_inf_proxy <- -1e300
max_inf_proxy <- 1e300

#file.dir <- dirname(parent.frame(2)$ofile)
setwd(dir)
getwd()


###################################################################################
## Load functions related to Udell's method 
###################################################################################
source("../Helper/helper_udell.R")

###################################################################################
## Load functions for data generation
###################################################################################
source("../Helper/helper_dataGeneration.R")
source("../Helper/helper_dataGeneration_linear.R")
source("../Helper/helper_dataGeneration_nonlinear.R")

###################################################################################
## Load functions for propensity score estimation and IPW and DR methods
###################################################################################
source("../Helper/helper_causalInference.R")
source("../Utils/miss.saem.v2.R")

###################################################################################
## Load functions for imputation
###################################################################################
source("../Helper/helper_imputation.R")


###################################################################################
## Load functions for simulations
###################################################################################
source("../Helper/helper_simulations.R")



###################################################################################
## Run DR simulations on incomplete data
###################################################################################
imp.methods <- c("mice", "mean", "mia.grf", "mean.grf", "mia.grf.ate", "mean.grf.ate", "saem", "udell")
results_compare_dr_miss <- data.frame()
results_compare_ipw_miss <- data.frame()
pp <- c(10)
r <- 3

nn <- c(100, 500, 1000, 5000)

for (cit in c(FALSE, TRUE)){
  for (cio in c(FALSE, TRUE)){
    for (p in pp){
      for (setting in c("linear1", "linear2")) {
        writeLines(paste0("Setting ", setting, "  CIT: ", cit, " CIO: ", cio, "\n"))
        for (ps.dependence in c("strong")) {
          for (use.mask in c(FALSE, TRUE)){
            for (n in nn) {
              writeLines(paste0("p: ", p,", n: ",n, ", use.mask: ", use.mask))
              results <- c()
              results$results_dr_miss <- c()
              results$results_ipw_miss <- c()
              try(results <- ate_estimation_miss(N=N, n=n, p = p, r = r,
                                                 prob = prob.miss, setting = setting, 
                                                 ps.dependence = ps.dependence,
                                                 nb.cores = nb_cores,
                                                 missing = missing,
                                                 trimming_weight = 1, 
                                                 use.mask = use.mask,
                                                 imputation.methods = imp.methods, 
                                                 cit = cit, cio = cio,
                                                 use.interaction = FALSE,
                                                 mi.m = mim))
              results_compare_dr_miss <- rbind(results_compare_dr_miss, 
                                               results$results_dr_miss)
              results_compare_ipw_miss <- rbind(results_compare_ipw_miss, 
                                                results$results_ipw_miss)
            }
          }
        }
      }
      save(results_compare_dr_miss,
           results_compare_ipw_miss,
           N, p, r, 
           alpha.star.strong, alpha.star.moderate, alpha.star.low, 
           tau, beta.star,
           prob.miss,
           file = paste0("./RData/", Date, "_N", N, "_p", p,"_nmin", min(nn), "_nmax", max(nn), "_", 
                         missing, prob.miss, "_results_compare",
                         "_cit", as.integer(cit), "_cio", as.integer(cio), "_miss.RData"))
    }
      
  }
  
}
