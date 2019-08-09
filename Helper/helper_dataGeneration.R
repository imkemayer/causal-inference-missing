generate_sim <- function(n, p = 10, r = 3, ng = 5, 
                         sd=0.1, 
                         seed = 0,
                         mechanism=FALSE, prop.missing = 0, 
                         cit = FALSE, cio = FALSE,
                         setting = "latentclass",
                         class.interaction = FALSE,
                         ps.dependence = "strong",
                         link = "nonlinear",
                         sigma.structure = "diagonal",
                         V = NULL){
  
  if (startsWith(setting, "linear")){
    sample <- gen_linear(n=n, p=p, r=r, setting = setting,
                         ps.dependence = ps.dependence,
                         sd=sd, seed=seed, 
                         mechanism=mechanism, prop.missing=prop.missing, 
                         cit=cit, cio=cio,
                         V = V)
    return(sample)
  }
  
  if (setting == "latentclass"){
    sample <- gen_latentclass(n, p, r, 
                              sd=sd, seed=seed, 
                              mechanism=mechanism, prop.missing=prop.missing, 
                              cit=cit, cio=cio,
                              class.interaction = class.interaction,
                              link = link)
    return(sample)
  }
  if (setting == "multisvd"){
    sample <- gen_multisvd(n, p, ng, 2, 2, 
                         sd=sd, seed=seed, 
                         mechanism=mechanism, prop.missing=prop.missing, 
                         cit=cit, cio=cio,
                         class.interaction = class.interaction,
                         link = link)
    return(sample)
  }
  if (setting == "dlvm"){
    sample <- gen_dlvm(n, p, d=3, 
                       sd=sd, seed=seed, 
                       mechanism=mechanism, prop.missing=prop.missing, 
                       cit=cit, cio=cio,
                       link = link,
                       sigma.structure = sigma.structure)
    return(sample)
  }
  
  if (setting == "dlvm2"){
    sample <- gen_dlvm2(n, p, d=3, 
                       sd=sd, seed=seed, 
                       mechanism=mechanism, prop.missing=prop.missing, 
                       cit=cit, cio=cio,
                       link = link,
                       sigma.structure = sigma.structure)
    return(sample)
  }
  
  if (setting == "ding1"){
    sample <- gen_ding(n, set = 1, 
                       seed = seed,
                       missing = "MNAR")
    return(sample)
  }
}  

expit <- function(x){
  return(1/(1+exp(-x)))
}
