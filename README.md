# causal-inference-missing

This repository contains codes and pipelines associated with the article [Doubly robust treatment effect estimation with missing attributes](http://dx.doi.org/10.1214/20-AOAS1356) by Mayer et al. (2020).

## General use
A full pipeline for estimating treatment effects in the presence of missing attributes, i.e., incomplete confounders and covariates, is provided in [`pipeline_causal_inference_with_missing_attributes.Rmd`](https://github.com/imkemayer/causal-inference-missing/blob/master/pipeline_causal_inference_with_missing_attributes.Rmd).
This pipeline can be applied directly on a custom data set (the default is a simulated toy example), provided that it suits the format as follows:

- `X.na`: confounders. A data.frame of size `#observations x #covariates`. With or without missing values.
- `W`: treatment assignment. A binary vector coded either with `{0,1}` or with `{FALSE,TRUE}` (representing `{control,treatment}`). Without missing values.
- `Y`: observed outcome. A numerical or binary vector (if binary, then coded with `{0,1}`). Without missing values.


## Application: Effect of Tranexamic Acid on Traumatic Brain Injury
The methodology has been applied on a medical question, the effect of the drug tranexamic acid on mortality among traumatic brain injury patients. The data used for this application is extracted from the [Traumabase® registry](http://www.traumabase.eu/en_US). This registry is only available upon request. However we provide the code used to analyse the data and to estimate the ATE in this context in [`TranexamicAcid/ate_analysis_traumabase_example.Rmd`](https://github.com/imkemayer/causal-inference-missing/blob/master/TranexamicAcid/ate_analysis_traumabase_example.Rmd).


## Reference:
Mayer, Imke, Erik Sverdrup, Tobias Gauss, Jean-Denis Moyer, Stefan Wager, and Julie Josse. 2020. "Doubly Robust Treatment Effect Estimation with Missing Attributes." Annals of Applied Statistics 14 (3): 1409–31.
