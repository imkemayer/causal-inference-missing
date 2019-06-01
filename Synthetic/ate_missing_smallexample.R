################################################################
# Example script to run a simple experiment with synthetic data 
# to compare different ATE estimators that handle missing values
################################################################

script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)
source("../Helper/helper_simulations.R")

# Proportion of missing values
prop.missing <- 0.5

# Choose methods to compare
imp.methods <- c("mice", "mean", "mia.grf", "mean.grf", "mia.grf.ate", "mean.grf.ate", "saem", "udell")

# Number of imputations for multiple imputations
mi.m <- 10

# Choose number of covariates
p <- 10

# For latent confounders settings, choose number of confounders (r <= p)
r <- 3

# Sample sizes
nn <- c(100, 500)

# Number of experiments per setting
N <- 2


results_compare_dr_miss <- data.frame()
results_compare_ipw_miss <- data.frame()

for (cit in c(FALSE, TRUE)){
  for (cio in c(FALSE, TRUE)){
  	for (use.mask in c(FALSE, TRUE)){
      for (n in nn) {
        writeLines(paste0("p: ", p,", n: ",n, ", use.mask: ", use.mask))
        results <- c()
        results$results_dr_miss <- c()
        results$results_ipw_miss <- c()
        try(results <- ate_estimation_miss(N=N, n=n, p = p, r = r,
                                           prob = prop.missing, setting = "latentclass", 
                                           missing = "MCAR", 
                                           use.mask = use.mask,
                                           imputation.methods = imp.methods, 
                                           cit = cit, cio = cio,
                                           mi.m = mi.m))
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
     prop.missing,
     file = paste0("./RData/", Sys.Date(), "_N", N, "_p", p,"_nmin", min(nn), "_nmax", max(nn), "_propNA", 
                    prop.missing, "_results_compare_miss.RData"))

