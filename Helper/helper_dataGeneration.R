generate_sim <- function(n, p=10, r=3, ng=5, 
                         tau=1,
                         sd=0.1, 
                         seed=0,
                         mechanism=FALSE, prop.missing=0, 
                         cit=FALSE, cio=FALSE,
                         cit2=FALSE, cio2=FALSE,
                         ci2_imp="mice",
                         setting="latentclass",
                         class.interaction=FALSE,
                         ps.dependence="strong",
                         link="nonlinear",
                         sigma.structure="diagonal",
                         V=NULL){
  
  if (startsWith(setting, "linear")){
    sample <- gen_linear(n=n, p=p, r=r, setting=setting,
                         tau=tau,
                         ps.dependence=ps.dependence,
                         sd=sd, seed=seed, 
                         mechanism=mechanism, prop.missing=prop.missing, 
                         cit=cit, cio=cio,
                         cit2=cit2, cio2=cio2,
                         V=V)
    return(sample)
  }
  
  if (setting == "latentclass"){
    sample <- gen_latentclass(n, p, r, 
                              sd=sd, seed=seed, 
                              mechanism=mechanism, prop.missing=prop.missing, 
                              cit=cit, cio=cio,
                              class.interaction=class.interaction,
                              link=link)
    return(sample)
  }
  if (setting == "multisvd"){
    sample <- gen_multisvd(n, p, ng, 2, 2, 
                         sd=sd, seed=seed, 
                         mechanism=mechanism, prop.missing=prop.missing, 
                         cit=cit, cio=cio,
                         class.interaction=class.interaction,
                         link=link)
    return(sample)
  }
  if (setting == "dlvm"){
    sample <- gen_dlvm(n, p, d=3, 
                       sd=sd, seed=seed, 
                       mechanism=mechanism, prop.missing=prop.missing, 
                       cit=cit, cio=cio,
                       link=link,
                       sigma.structure=sigma.structure,
                       cit2=cit2, cio2=cio2,
                       ci2_imp=ci2_imp)
    return(sample)
  }
  
  if (setting == "dlvm2"){
    sample <- gen_dlvm2(n, p, d=3, 
                       sd=sd, seed=seed, 
                       mechanism=mechanism, prop.missing=prop.missing, 
                       cit=cit, cio=cio,
                       link=link,
                       sigma.structure=sigma.structure)
    return(sample)
  }
  
  if (setting == "ding1"){
    sample <- gen_ding(n, set=1,tau=tau,
                       seed=seed,
                       missing="MNAR")
    return(sample)
  }
}  

expit <- function(x){
  return(1/(1+exp(-x)))
}
