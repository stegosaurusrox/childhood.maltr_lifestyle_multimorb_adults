# CREDIT:
# some of this code was written by
# Serena Defina: https://github.com/SereDef
# and adjusted, merged, extended by Esther Walton

## Prep the environment 

rm(list=ls())

# Load needed packages
pack <- c('mice','stringr','miceadds','nnet','MNLpred','openxlsx')
invisible(lapply(pack, require, character.only = T));


# Date (for marking output file)
date <- format(Sys.Date(), "%d%m%y")
mod_imp <- readRDS(file = ".\\data\\adult.moderator_data_imp.rds")

# ------------------------------------------------------------------------------

# name cohort
cohort <- 'ALSPAC' # UKBB, ALSPAC


# Add standardized scores for main variables
mods_partners = c('PhyAct','Smoke','Alcohol', 'Diet_healthy', 'Diet_trad', 'Diet_proc.conf', 'Diet_semiveget')
mods_mothers = c('PhyAct','Smoke','Alcohol', 'Diet_healthy', 'Diet_proc', 'Diet_conf', 'Diet_veget')
mods_UKBB = c('PhyAct','Smoke','Alcohol')

mods=mods_partners

# define covariates (with or without sex; and then either PCM.fu without baseline correction - or PCM with, but not both)
covs_wSex <- '+ Age + Educ + Ethn + Sex'
covs_wSex_fu <- '+ Age.fu + Educ + Ethn + Sex'
covs_ALSPAC <- '+ Age + Educ + Ethn'
covs_ALSPAC_baseline.MDD <- '+ Age + Educ + Ethn + MDD.baseline'

covs=covs_ALSPAC_baseline.MDD

#====
  
to_trans <- c(mods)
mod_imp <- datlist2mids( scale_datlist( mids2datlist(mod_imp), 
                             orig_var = to_trans, trafo_var = paste0(to_trans, '_z')))

# Create numeric CM variables (for plotting and pps later)
longdat <- complete(mod_imp,'long', include=T)
longdat$CM_num <- as.numeric(levels(longdat$CM))[longdat$CM]
mod_imp <- as.mids(longdat)

# Save finished mids 
#saveRDS(mod_imp, file.path(datapath,'imp_sample.rds'))

# Save each imputed setfor portenial plottining
#for (m in 0) { write.csv(complete(mod_imp,m), file = paste0('./data/byimp/imp',m,'.csv')) }

# ------------------------------------------------------------------------------


pool_fit <- function(mods=mods,
                     z_or_bin='_z', interact='*',
                     outc='PCM', exp='CM') {
  cat(outc, '\n')
  # Define covariates
  covs <- covs
  # Initiate stack
  models = data.frame()
  
  for (modr in paste0(mods,z_or_bin)) {
     
      fit <- with(mod_imp, nnet::multinom(as.formula(paste(outc,'~',exp, interact, modr, covs)), 
                                              model=T, trace=F));
    p_fit <- mice::pool(fit) # pool results 
    mod <- summary(p_fit) # extract relevant information
    mod[,-c(1,2)] <- round(mod[,-c(1,2)],4)
    mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
    
    #cat('--------', as.character(mod[nrow(mod),'term']), '--> p =', round(mod$p.value[nrow(mod)], 3), '\n')
    
    
      if (endsWith(outc, 'REC')) { levels(mod$y.level) <- c("C:healthy","C:intern", "C:adipos")
      } else { levels(mod$y.level) <- c("H:intern", "H:adipos", "H:comorb") }
      mod$OR  <- round(exp(mod$estimate), 4)
      mod$lci <- round(exp((mod$estimate) - 1.96*mod$std.error), 4)
      mod$uci <- round(exp((mod$estimate) + 1.96*mod$std.error), 4)
      mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
    
    mod <- cbind(rep(modr, nrow(mod)),mod)
    
    models <- rbind(models, mod, rep(NA, ncol(mod)))
  }
  names(models)[1] = 'model' # print(mod)
  return(models)
}
# ------------------------------------------------------------------------------

# Baseline analyses (comorbidity [vs. healthy], 

CM_modz_PCM_int <- pool_fit(exp = 'CM',
                    mods = mods,
                    z_or_bin = '_z',
                    interact = '*',
                    outc = 'PCM')


CM_modb_PCM_int <- pool_fit(exp = 'CM',
                            mods = mods,
                            z_or_bin = '_bin',
                            interact = '*',
                            outc = 'PCM')


CM_modz_PCM_noint <- pool_fit(exp = 'CM',
                            mods = mods,
                            z_or_bin = '_z',
                            interact = '+',
                            outc = 'PCM')


CM_modb_PCM_noint <- pool_fit(exp = 'CM',
                            mods = mods,
                            z_or_bin = '_bin',
                            interact = '+',
                            outc = 'PCM')


# follow-up (ALSPAC partners, mothers and UKBB)
CM_modz_PCM_int <- pool_fit(exp = 'CM',
                            mods = mods,
                            z_or_bin = '_z',
                            interact = '*',
                            outc = 'PCM.fu')


CM_modb_PCM_int <- pool_fit(exp = 'CM',
                            mods = mods,
                            z_or_bin = '_bin',
                            interact = '*',
                            outc = 'PCM.fu')


CM_modz_PCM_noint <- pool_fit(exp = 'CM',
                              mods = mods,
                              z_or_bin = '_z',
                              interact = '+',
                              outc = 'PCM.fu')


CM_modb_PCM_noint <- pool_fit(exp = 'CM',
                              mods = mods,
                              z_or_bin = '_bin',
                              interact = '+',
                              outc = 'PCM.fu')

# ------------------------------------------------------------------------------
names = ls()[grepl('PCM', ls())]
for (n in names) { cat("'",n,"'= ",n,", ", sep='')}

modls <- list('CM_modb_PCM_int'= CM_modb_PCM_int, 'CM_modb_PCM_noint'= CM_modb_PCM_noint,
              'CM_modz_PCM_int'= CM_modz_PCM_int, 'CM_modz_PCM_noint'= CM_modz_PCM_noint)

openxlsx::write.xlsx(modls, file = paste0('./results/',cohort,'_Results_imp_',date,'.xlsx'), 
                     overwrite=T)

