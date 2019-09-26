.libPaths( c( .libPaths(), "~/R/R_libs/") )
previous.wd <- getwd()
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)
try(source("../Helper/helper_dataGeneration.R"))
try(source("../Helper/helper_dataGeneration_linear.R"))
try(source("../Helper/helper_dataGeneration_nonlinear.R"))
try(source("../Helper/helper_causalInference.R"))
try(source("../Helper/helper_imputation.R"))
try(source("../Helper/helper_udell.R"))
try(source("../Helper/helper_simulations.R"))
try(source("./Helper/helper_dataGeneration.R"))
try(source("./Helper/helper_dataGeneration_linear.R"))
try(source("./Helper/helper_dataGeneration_nonlinear.R"))
try(source("./Helper/helper_causalInference.R"))
try(source("./Helper/helper_imputation.R"))
try(source("./Helper/helper_udell.R"))
try(source("./Helper/helper_simulations.R"))
try(source("./helper_dataGeneration.R"))
try(source("./helper_dataGeneration_linear.R"))
try(source("./helper_dataGeneration_nonlinear.R"))
try(source("./helper_causalInference.R"))
try(source("./helper_imputation.R"))
try(source("./helper_udell.R"))
try(source("./helper_simulations.R"))
try(source("../Utils/miss.saem.v2.R"))
try(source("./Utils/miss.saem.v2.R"))
try(source("../Utils/amputation.R"))
try(source("./Utils/amputation.R"))

setwd(previous.wd)

if (!require(doParallel)) {
	install.packages('doParallel')
	library(doParallel)
}
# if (!require(doSNOW)) {
#   install.packages('doSNOW')
#   library(doSNOW)
# }
if (!require(foreach)) {
	install.packages('foreach')
	library(foreach)
}

if (!require(missMDA)) {
  install.packages('missMDA')
  library(missMDA)
}
if (!require(norm)) {
  install.packages('norm')
  library(norm)
}
if (!require(MASS)) {
  install.packages('MASS')
  library(MASS)
}
if (!require(misaem)) {
  install.packages('misaem')
  library(misaem)
}
if (!require(grf)) {
  install.packages('grf')
  library(grf)
}
if (!require(caret)) {
  install.packages('caret')
  library(caret)
}
if (!require(ranger)) {
  install.packages('ranger')
  library(ranger)
}
if (!require(FactoMineR)) {
  install.packages('FactoMineR')
  library(FactoMineR)
}
if (!require(dplyr)) {
  install.packages('dplyr')
  library(dplyr)
}

if (!require(assertthat)) {
	install.packages('assertthat')
	library(assertthat)
}
if (!require(pracma)) {
	install.packages('pracma')
	library(pracma)
}
if (!require(softImpute)) {
	install.packages('softImpute')
	library(softImpute)
}

if (!require(spcov)) install.packages('spcov')
if (!require(parallel)) install.packages('parallel')
if (!require(mice)) install.packages('mice')
if (!require(dplyr)) install.packages('dplyr')
# if (!require(dbarts)) install.packages('dbarts')

