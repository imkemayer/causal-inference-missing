# Numeric variables (will be rounded to 1 decimal digit)
data.num <- c("Lactates", "Fibrinogene.1", "Hb", "Alcool","Dose.NAD.depart", "DTC.IP.max")

# Imputation using FAMD
get_FAMD <- function(database,seed, 
                     ncp=5, Method = "Regularized", scale = T, threshold = 1e-06) {
  imputed <- imputeFAMD(database, ncp=ncp,seed = seed)
  imputedData <- imputed$completeObs
  
  for (j in colnames(imputedData)){
    imputedData[,j] <- cast_types(j,imputedData)
  }
  
  # Re-impose assumption that no mydriase at the pre-hospital phase implies that there was no mannitol given:
  if ("Mannitol.SSH" %in% colnames(imputedData)){
    imputedData[imputedData$Mydriase %in% c("Non"),"Mannitol.SSH"] <- "No mydriase"
  }
  
  # Re-imposing the "Not tested" in Regression.mydriase.sous.osmotherapie
  if ("Regression.mydriase.sous.osmotherapie" %in% colnames(imputedData)){
    imputedData[imputedData$Mydriase %in% c("Non"),"Regression.mydriase.sous.osmotherapie"] <- "Not tested"
    imputedData[imputedData$Mannitol.SSH %in% c("Rien"),"Regression.mydriase.sous.osmotherapie"] <- "Not tested"
  }
  
  return(imputedData)
}

# Recode values of imputed categorical variables and recast some numericals variables into integers 
cast_types = function(i,df){
  if (is.factor(df[,i])){
    df[,i] = plyr::mapvalues(df[,i], from = levels(df[,i]), to = gsub(paste(i,"_",sep=''), "", levels(df[,i])))
  } else {
    if(i %in% data.num){
      df[,i] <- round(df[,i],digits=1)
    }  else{
      df[,i] <- as.integer(round(df[,i],digits=0))
    }
  }
  return(df[,i])
}




################################################################################################################

# Using SAEM for PS estimation we need to recode qualitative variables into numerical.
prepare_data_misaem <- function(df, use_targets=T){
  recode_covariates_1 <- c("Trauma.Center","Traitement.anticoagulant", "Traitement.antiagregants",
                           "Glasgow.initial", "Glasgow.moteur.initial" ,
                           "ACR.1", "Hemocue.init", "Delta.hemocue", "Shock.index.ph", "Delta.shock.index",
                           "Catecholamines", 
                           "IOT.SMUR", 
                           "Ventilation.FiO2", 
                           "KTV.poses.avant.TDM", 
                           "Choc.hemorragique","Trauma.cranien",
                           "PIC",
                           "DVE", "Hypothermie.therapeutique", "Craniectomie.decompressive",
                           "Acide.tranexamique", "Bloc.J0.neurochirurgie",
                           "DC.Trauma", 
                           "Regression.mydriase.sous.osmotherapie")
  
  # if ("Regression.mydriase.sous.osmotherapie" %in% colnames(df)){
  #   recoded_df <- df %>%
  #                   dplyr::select(-"Regression.mydriase.sous.osmotherapie")
  # } else {
  #   recoded_df <- df
  # }
  
  recoded_df <- df
  if ("Regression.mydriase.sous.osmotherapie" %in% colnames(df)){
    levels(recoded_df$Regression.mydriase.sous.osmotherapie) <- c("0", "1", "-1")
  }
  
  if ("Trauma.Center" %in% colnames(df)){
    levels(recoded_df$Trauma.Center) <- 1:length(levels(recoded_df$Trauma.Center))
  }
  
  for (j in recode_covariates_1){
    if (j %in% colnames(recoded_df)){
      recoded_df[,j] <- as.numeric(as.character(recoded_df[,j]))
    }
  }
  
  
  # Recode osmotherapy treatment into binary (and then numeric)
  
  if ("Mannitol.SSH" %in% colnames(recoded_df)){
    new_Mannitol.SSH <- rep(0, length(recoded_df$Mannitol.SSH))
    new_Mannitol.SSH[recoded_df$Mannitol.SSH %in% c("Mannitol", "SSH")] <- 1
    new_Mannitol.SSH[is.na(recoded_df$Mannitol.SSH)] <- NA
    recoded_df$Mannitol.SSH <- new_Mannitol.SSH
  }
  
  if ("Osmotherapie" %in% colnames(recoded_df)){
    new_Osmotherapie <- rep(0, length(recoded_df$Osmotherapie))
    new_Osmotherapie[recoded_df$Osmotherapie %in% c("Mannitol", "SSH")] <- 1
    new_Osmotherapie[is.na(recoded_df$Osmotherapie)] <- NA
    recoded_df$Osmotherapie <- new_Osmotherapie
  }
  
  # Recode mydriase/anomalie pupillaire into binary (and then numeric)
  if ("Mydriase" %in% colnames(recoded_df)){
    new_Mydriase <- rep(0, length(recoded_df$Mydriase))
    new_Mydriase[recoded_df$Mydriase %in% c("Anisocorie (unilatérale)", "Anomalie pupillaire")] <- 1
    new_Mydriase[is.na(recoded_df$Mydriase)] <- NA
    recoded_df$Mydriase <- new_Mydriase
  }
  
  if ("Anomalie.pupillaire" %in% colnames(recoded_df)){
    new_Anomalie.pupillaire <- rep(0, length(recoded_df$Anomalie.pupillaire))
    new_Anomalie.pupillaire[recoded_df$Anomalie.pupillaire %in% c("Anisocorie (unilatérale)", "Mydriase (bilatérale)")] <- 1
    new_Anomalie.pupillaire[is.na(recoded_df$Anomalie.pupillaire)] <- NA
    recoded_df$Anomalie.pupillaire <- new_Anomalie.pupillaire
  }
  
  
  if ("Glasgow.sortie" %in% colnames(recoded_df)) {
    recoded_df <- recoded_df %>%
      dplyr::select(-"Glasgow.sortie")
  }
  
  if (use_targets){
    recoded_df <- recoded_df %>%
      filter(Trauma.cranien == 1 | AIS.tete >= 2)
  }
  
  return(recoded_df)
}
