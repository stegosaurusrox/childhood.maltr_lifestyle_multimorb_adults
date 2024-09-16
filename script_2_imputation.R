# CREDIT:
# some of this code was written by
# Camille Souama: https://github.com/camillesouama/earlycause-tools/tree/main/Amsterdam%20UMC/Meta-analysis%20on%20childhood%20maltreatment%20and%20(comorbid)%20depression%20and%20cardiometabolic%20disease
# Serena Defina: https://github.com/SereDef
# and adjusted, merged, extended by Esther Walton

library(dplyr)
library(mice)


# Date (for marking output file)
date <- format(Sys.Date(), '%d%m%y')

# ------------------------------------------------------------------------------
# read in file (unless already loaded)
Data <- readRDS(file = './data/adult.moderator_data_raw.rds')

# for UKBB: delete unneeded variables
var=c('CM_sev2','CIDI.MDD.No.Info','CIDI.MDD.Screen','CIDI.MDD.Response',
      'f.20514.0.0','f.20507.0.0','f.20510.0.0','f.20508.0.0','f.20517.0.0','f.20518.0.0','f.20519.0.0','f.20513.0.0','f.20511.0.0',
      'PHQ9.Severity',
      'V1', 'SRDepression' ,'InterviewDepression' ,'DrugDepression' ,#'AntidepMed' ,
      'InterviewDiabetesType1' ,'InterviewDiabetesType2',
      'DiabMed_women_baseline', 'CVDMed_women_baseline', 'DiabMed_men_baseline', 'CVDMed_men_baseline', 'Depressed.Ever',
      'PhyAct_z', 'Smoke_z','Alcohol_z')

Data$DM.fu=as.factor(Data$DM.fu)

Data = Data[,!(names(Data) %in% var)]

## for all: delete all higher level variables before imputing lower-level vars

# ALSPAC partners
higher.vars=c('CM','CM_sev',
              'MDD.baseline','MDD.fu',
              'DM.fu','CMD','CMD.fu',
              'PCM','PCM.fu',
              'PhyAct_bin','Smoke_bin','Alcohol_bin',
              'Diet_healthy_bin','Diet_trad_bin','Diet_proc.conf_bin','Diet_semiveget_bin')

# ALSPAC mothers
higher.vars=c('qlet','mult',
              'CM','CM_sev',
              'MDD.baseline','MDD.fu',
              'CMD','CardiacMD','CVD2MD','CVD2',
              'CD', # tried in mothers to also remove CD, as this might be driving the very high CVD1 rates
              'PCM','PCM.fu',
              'Smoke_bin','Alcohol_bin',
              'Diet_healthy_bin','Diet_conf_bin','Diet_proc_bin','Diet_veget_bin')

# UKBB
higher.vars=c('CM','CM_sev',
              'MDD_med',
              'PM.fu','PM',
              'PCV.fu','PCV',
              'PCM_CVD2.fu','PCM_CVD2','PCM_CD.fu','PCM_CD',
              'CVD1.fu','CVD2.fu',
              'CMD','CMD.fu',
              'CardiacMD','CardiacMD.fu',
              'CVD2MD','CVD2MD.fu',
              'CD.fu', 
              'PCM','PCM.fu',
              'Alcohol_bin')

###########################
Data = Data[,!(names(Data) %in% higher.vars)]

# ==============================================================================
# IMPUTATION ===================================================================
# ==============================================================================
# Initiate


# Make a mice method object 
# NOTE: use random forest for continuous data, logistic regression for binomial (sex) and polyreg for categorical
meth <- mice::make.method(Data, defaultMethod = c("rf", "logreg", "polyreg"))

# Random forest imputation (ran in parallel)
mod_imp <- mice::futuremice(Data, method = meth, m = 20, maxit = 40, ntree = 10, rfPackage = 'ranger',
                           parallelseed = 310896,
                           print=TRUE)



# add back in higher-level variables



sample= 'ALSPAC_mothers' # 'ALSPAC_partners', 'ALSPAC_mothers' 'UKBB'

if(sample=='ALSPAC_partners'){
  
  ################# ALSPAC partners
  
 
  longdat <- complete(mod_imp,'long', include=T)
  
  # Code childhood maltreatment variables : 1=yes, 0=no
  

  
  longdat$CM <- as.factor(ifelse(longdat$CM_PA == 1 | longdat$CM_EA == 1 | longdat$CM_SA == 1, 1,
                                 ifelse(longdat$CM_PA == 0 & longdat$CM_EA == 0 & longdat$CM_SA == 0, 0, NA)))
  
  longdat$CM_sev <- as.factor(rowSums(data.frame(as.numeric(as.character(longdat$CM_PA)),
                                                 as.numeric(as.character(longdat$CM_EA)),
                                                 as.numeric(as.character(longdat$CM_SA)))))
  
  # Code depression-related variables : 1=yes, 0=no 
  
  # MDD at any baseline
  var_MDD=c('MDD.baseline1', 'MDD.baseline2', 'MDD.baseline3', 'MDD.baseline4')
  longdat$MDD.baseline <- ifelse(rowSums(longdat[,var_MDD], na.rm = T)>0 , 1, 
                                 ifelse(rowSums(is.na(longdat[,var_MDD]))==length(var_MDD), NA, 0))
  
  # MDD at follow-up
  longdat$MDD.fu <- ifelse(longdat$MDD.baseline==0 & longdat$MDD.fu1==1, 1, 
                           ifelse(rowSums(longdat[,c(var_MDD,'MDD.fu1')], na.rm = T)==0, 0, NA))
  
  longdat$MDD.fu=as.factor(longdat$MDD.fu)
  longdat$MDD.fu=as.factor(longdat$MDD.fu)
  longdat$MDD.fu1 <- as.factor(longdat$MDD.fu1)
  
  longdat$DM.fu <- ifelse(longdat$DM.baseline==0 & longdat$DM.fu1==1, 1, 
                          ifelse(longdat$DM.fu1==0, 0, NA))
  
  # Code cardiometabolic-related variables: 1=yes, 0=no
  
  longdat$DM.fu=as.factor(longdat$DM.fu)
  
  # Code PCM multimorbidity variables: 0=healthy, 1=MDD, 2=CMD, 3=comorbid
  longdat$CMD <- as.factor(ifelse(longdat$CVD1 == 0 & longdat$DM.fu1 == 0,0,
                                  ifelse(longdat$CVD1 == 1 | longdat$DM.fu1 == 1,1,NA)))
  
  longdat$CMD.fu <- as.factor(ifelse(longdat$CVD1 == 0 & longdat$DM.fu == 0,0,
                                     ifelse(longdat$CVD1 == 1 | longdat$DM.fu == 1,1,NA)))
  
  longdat$PCM <- as.factor(ifelse(longdat$MDD.fu1 == 0 & longdat$CMD == 0, 0,
                                  ifelse(longdat$MDD.fu1 == 1 & longdat$CMD == 0, 1,
                                         ifelse(longdat$MDD.fu1 == 0 & longdat$CMD == 1, 2,
                                                ifelse(longdat$MDD.fu1 == 1 & longdat$CMD == 1, 3, NA)))))
  
  
  longdat$PCM.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$CMD == 0, 0,
                                     ifelse(longdat$MDD.fu == 1 & longdat$CMD == 0, 1,
                                            ifelse(longdat$MDD.fu == 0 & longdat$CMD == 1, 2,
                                                   ifelse(longdat$MDD.fu == 1 & longdat$CMD == 1, 3, NA)))))
  
  # binary exercise following WHO guidelines of >3 times p/w 
  
  longdat$PhyAct_bin = as.factor(ifelse(longdat$PhyAct>4,1,
                                        ifelse(longdat$PhyAct<5,0,NA)))
  
  
  # binary smoking
  longdat$Smoke_bin <- as.factor(ifelse(longdat$Smoke>0,1,
                                        ifelse(longdat$Smoke==0,0,NA)))
  
  
  # binary alcohol
  longdat$Alcohol_bin <- as.factor(ifelse(longdat$Alcohol>4,1,
                                          ifelse(longdat$Alcohol<5,0,NA)))
  
  
  # binary diet
  # median split
  longdat$Diet_healthy_bin=as.factor(ifelse(longdat$Diet_healthy>median(longdat$Diet_healthy, na.rm=T),1,
                                            ifelse(longdat$Diet_healthy<=median(longdat$Diet_healthy, na.rm=T),0,NA)))
  
  longdat$Diet_trad_bin=as.factor(ifelse(longdat$Diet_trad>median(longdat$Diet_trad, na.rm=T),1,
                                         ifelse(longdat$Diet_trad<=median(longdat$Diet_trad, na.rm=T),0,NA)))
  
  longdat$Diet_proc.conf_bin=as.factor(ifelse(longdat$Diet_proc.conf>median(longdat$Diet_proc.conf, na.rm=T),1,
                                              ifelse(longdat$Diet_proc.conf<=median(longdat$Diet_proc.conf, na.rm=T),0,NA)))
  
  longdat$Diet_semiveget_bin=as.factor(ifelse(longdat$Diet_semiveget>median(longdat$Diet_semiveget, na.rm=T),1,
                                              ifelse(longdat$Diet_semiveget<=median(longdat$Diet_semiveget, na.rm=T),0,NA)))
  
} else if (sample=='ALSPAC_mothers'){

  #################################### ALSPAC mothers
  
  

  longdat <- complete(mod_imp,'long', include=T)
  
  # Code childhood maltreatment variables : 1=yes, 0=no
  
   longdat$CM <- as.factor(ifelse(longdat$CM_PA == 1 | longdat$CM_EA == 1 | longdat$CM_SA == 1, 1,
                                 ifelse(longdat$CM_PA == 0 & longdat$CM_EA == 0 & longdat$CM_SA == 0, 0, NA)))
  
  longdat$CM_sev <- as.factor(rowSums(data.frame(as.numeric(as.character(longdat$CM_PA)),
                                                 as.numeric(as.character(longdat$CM_EA)),
                                                 as.numeric(as.character(longdat$CM_SA)))))
  
  # Code depression-related variables : 1=yes, 0=no 
  
  # MDD at any baseline
  var_MDD=c('MDD.baseline1', 'MDD.baseline2', 'MDD.baseline3', 'MDD.baseline4', 'MDD.baseline5')
  longdat$MDD.baseline <- ifelse(rowSums(longdat[,var_MDD], na.rm = T)>0 , 1, 
                                 ifelse(rowSums(is.na(longdat[,var_MDD]))==length(var_MDD), NA, 0))
  
  # MDD at follow-up
  longdat$MDD.fu <- ifelse(longdat$MDD.baseline==0 & longdat$MDD.fu1==1, 1, 
                           ifelse(rowSums(longdat[,c(var_MDD,'MDD.fu1')], na.rm = T)==0, 0, NA))
  
  longdat$MDD.fu=as.factor(longdat$MDD.fu)
  longdat$MDD.fu1 <- as.factor(longdat$MDD.fu1)
  
  
  
  # Code PCM multimorbidity variables : 0=healthy, 1=MDD, 2=CMD, 3=comorbid
  longdat$CMD <- as.factor(ifelse(longdat$CVD1 == 0 & longdat$DM == 0,0,
                                  ifelse(longdat$CVD1 == 1 | longdat$DM == 1,1,NA)))
  
  
  longdat$PCM <- as.factor(ifelse(longdat$MDD.fu1 == 0 & longdat$CMD == 0, 0,
                                  ifelse(longdat$MDD.fu1 == 1 & longdat$CMD == 0, 1,
                                         ifelse(longdat$MDD.fu1 == 0 & longdat$CMD == 1, 2,
                                                ifelse(longdat$MDD.fu1 == 1 & longdat$CMD == 1, 3, NA)))))
  
  
  longdat$PCM.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$CMD == 0, 0,
                                     ifelse(longdat$MDD.fu == 1 & longdat$CMD == 0, 1,
                                            ifelse(longdat$MDD.fu == 0 & longdat$CMD == 1, 2,
                                                   ifelse(longdat$MDD.fu == 1 & longdat$CMD == 1, 3, NA)))))
  
  
  
  # 
  # 
  # # binary smoking
  longdat$Smoke_bin <- as.factor(ifelse(longdat$Smoke>0,1,
                                        ifelse(longdat$Smoke==0,0,NA)))
  
  
  # binary alcohol
  longdat$Alcohol_bin <- as.factor(ifelse(longdat$Alcohol>4,1,
                                          ifelse(longdat$Alcohol<5,0,NA)))
  
  
  # # binary diet
  # # median split
  longdat$Diet_healthy_bin=as.factor(ifelse(longdat$Diet_healthy>median(longdat$Diet_healthy, na.rm=T),1,
                                            ifelse(longdat$Diet_healthy<=median(longdat$Diet_healthy, na.rm=T),0,NA)))
  
  longdat$Diet_proc_bin=as.factor(ifelse(longdat$Diet_proc>median(longdat$Diet_proc, na.rm=T),1,
                                         ifelse(longdat$Diet_proc<=median(longdat$Diet_proc, na.rm=T),0,NA)))
  
  longdat$Diet_conf_bin=as.factor(ifelse(longdat$Diet_conf>median(longdat$Diet_conf, na.rm=T),1,
                                         ifelse(longdat$Diet_conf<=median(longdat$Diet_conf, na.rm=T),0,NA)))
  
  longdat$Diet_veget_bin=as.factor(ifelse(longdat$Diet_veget>median(longdat$Diet_veget, na.rm=T),1,
                                          ifelse(longdat$Diet_veget<=median(longdat$Diet_veget, na.rm=T),0,NA)))
  
} else if (sample==UKBB) {
  
  longdat <- complete(mod_imp,'long', include=T)
  
  longdat$CM <- as.factor(ifelse(longdat$CM_PA == 1 | longdat$CM_EA == 1 | longdat$CM_SA == 1, 1,
                              ifelse(longdat$CM_PA == 0 & longdat$CM_EA == 0 & longdat$CM_SA == 0, 0, NA)))
  
  
  longdat$CM_sev <- as.factor(rowSums(data.frame(as.numeric(as.character(longdat$CM_PA)),
                                  as.numeric(as.character(longdat$CM_EA)),
                                  as.numeric(as.character(longdat$CM_SA)))))

  longdat$MDD_med <- as.factor(ifelse((!is.na(longdat$MDD) & longdat$MDD == 0) & (is.na(longdat$AntidepMed)), 0,
                                   ifelse((!is.na(longdat$MDD) & longdat$MDD == 1) | (!is.na(longdat$AntidepMed) & longdat$AntidepMed == 1), 1, NA)))
  
  
  longdat <- longdat %>%
    mutate(
      CD.fu = case_when(
        CD == 0 & (CD.fu1 == 1 | CD.fu2 == 1) ~ 1,
        CD == 0 | CD.fu1 == 0 | CD.fu2 == 0 ~ 0
      )
    )
  
  longdat$CD.fu=as.factor(longdat$CD.fu)
  
  
  # CVD1.fu
  longdat <- longdat %>%
    mutate(
      CVD1.fu = case_when(
        CVD1 == 0 & (CVD1.fu1 == 1 | CVD1.fu2 == 1 | CVD1.fu3 == 1) ~ 1,
        CVD1 == 0 | CVD1.fu1 == 0 | CVD1.fu2 == 0 | CVD1.fu3 == 0 ~ 0
      )
    )
  
  
  longdat$CVD1.fu=as.factor(longdat$CVD1.fu)
  
  # CVD2.fu
  longdat <- longdat %>%
    mutate(
      CVD2.fu = case_when(
        CVD2 == 0 & (CVD2.fu1 == 1 | CVD2.fu2 == 1 | CVD2.fu3 == 1) ~ 1,
        CVD2 == 0 | CVD2.fu1 == 0 | CVD2.fu2 == 0 | CVD2.fu3 == 0 ~ 0
      )
    )
  
  
  longdat$CMD <- as.factor(ifelse(longdat$CVD1 == 0 & longdat$DM_baseline == 0,0,
                               ifelse(longdat$CVD1 == 1 | longdat$DM_baseline == 1,1,NA)))
  
  longdat$CMD.fu <- as.factor(ifelse(longdat$CVD1.fu == 0 & longdat$DM.fu == 0,0,
                                  ifelse(longdat$CVD1.fu == 1 | longdat$DM.fu == 1,1,NA)))
  
  longdat$CardiacMD <- as.factor(ifelse(longdat$CD == 0 & longdat$DM_baseline == 0,0,
                                     ifelse(longdat$CD == 1 | longdat$DM_baseline == 1,1,NA)))
  
  longdat$CardiacMD.fu <- as.factor(ifelse(longdat$CD.fu == 0 & longdat$DM.fu == 0,0,
                                        ifelse(longdat$CD.fu == 1 | longdat$DM.fu == 1,1,NA)))
  
  longdat$CVD2MD <- as.factor(ifelse(longdat$CVD2 == 0 & longdat$DM_baseline == 0,0,
                                  ifelse(longdat$CVD2 == 1 | longdat$DM_baseline == 1,1,NA)))
  
  longdat$CVD2MD.fu <- as.factor(ifelse(longdat$CVD2.fu == 0 & longdat$DM.fu == 0,0,
                                     ifelse(longdat$CVD2.fu == 1 | longdat$DM.fu == 1,1,NA)))
  
  
  longdat$PCM <- as.factor(ifelse(longdat$MDD == 0 & longdat$CMD == 0, 0,
                               ifelse(longdat$MDD == 1 & longdat$CMD == 0, 1,
                                      ifelse(longdat$MDD == 0 & longdat$CMD == 1, 2,
                                             ifelse(longdat$MDD == 1 & longdat$CMD == 1, 3, NA)))))
  
  longdat$PCM.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$CMD.fu == 0, 0,
                                  ifelse(longdat$MDD.fu == 1 & longdat$CMD.fu == 0, 1,
                                         ifelse(longdat$MDD.fu == 0 & longdat$CMD.fu == 1, 2,
                                                ifelse(longdat$MDD.fu == 1 & longdat$CMD.fu == 1, 3, NA)))))
  
  longdat$PCM_CD <- as.factor(ifelse(longdat$MDD == 0 & longdat$CardiacMD == 0, 0,
                                  ifelse(longdat$MDD == 1 & longdat$CardiacMD == 0, 1,
                                         ifelse(longdat$MDD == 0 & longdat$CardiacMD == 1, 2,
                                                ifelse(longdat$MDD == 1 & longdat$CardiacMD == 1, 3, NA)))))
  
  longdat$PCM_CD.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$CardiacMD.fu == 0, 0,
                                     ifelse(longdat$MDD.fu == 1 & longdat$CardiacMD.fu == 0, 1,
                                            ifelse(longdat$MDD.fu == 0 & longdat$CardiacMD.fu == 1, 2,
                                                   ifelse(longdat$MDD.fu == 1 & longdat$CardiacMD.fu == 1, 3, NA)))))
  
  longdat$PCM_CVD2 <- as.factor(ifelse(longdat$MDD == 0 & longdat$CVD2MD == 0, 0,
                                    ifelse(longdat$MDD == 1 & longdat$CVD2MD == 0, 1,
                                           ifelse(longdat$MDD == 0 & longdat$CVD2MD == 1, 2,
                                                  ifelse(longdat$MDD == 1 & longdat$CVD2MD == 1, 3, NA)))))
  
  longdat$PCM_CVD2.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$CVD2MD.fu == 0, 0,
                                       ifelse(longdat$MDD.fu == 1 & longdat$CVD2MD.fu == 0, 1,
                                              ifelse(longdat$MDD.fu == 0 & longdat$CVD2MD.fu == 1, 2,
                                                     ifelse(longdat$MDD.fu == 1 & longdat$CVD2MD.fu == 1, 3, NA)))))
  
  longdat$PCV <- as.factor(ifelse(longdat$MDD == 0 & longdat$CVD1 == 0, 0,
                               ifelse(longdat$MDD == 1 & longdat$CVD1 == 0, 1,
                                      ifelse(longdat$MDD == 0 & longdat$CVD1 == 1, 2,
                                             ifelse(longdat$MDD == 1 & longdat$CVD1 == 1, 3, NA)))))
  
  longdat$PCV.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$CVD1.fu == 0, 0,
                                  ifelse(longdat$MDD.fu == 1 & longdat$CVD1.fu == 0, 1,
                                         ifelse(longdat$MDD.fu == 0 & longdat$CVD1.fu == 1, 2,
                                                ifelse(longdat$MDD.fu == 1 & longdat$CVD1.fu == 1, 3, NA)))))
  
  longdat$PM <- as.factor(ifelse(longdat$MDD == 0 & longdat$DM_baseline == 0, 0,
                              ifelse(longdat$MDD == 1 & longdat$DM_baseline == 0, 1,
                                     ifelse(longdat$MDD == 0 & longdat$DM_baseline == 1, 2,
                                            ifelse(longdat$MDD == 1 & longdat$DM_baseline == 1, 3, NA)))))
  
  longdat$PM.fu <- as.factor(ifelse(longdat$MDD.fu == 0 & longdat$DM.fu == 0, 0,
                                 ifelse(longdat$MDD.fu == 1 & longdat$DM.fu == 0, 1,
                                        ifelse(longdat$MDD.fu == 0 & longdat$DM.fu == 1, 2,
                                               ifelse(longdat$MDD.fu == 1 & longdat$DM.fu == 1, 3, NA)))))
  
  
  # recoded to problematic drinking
  longdat$Alcohol_bin <- as.factor(ifelse(longdat$Alcohol <4, 0, 
                                       ifelse(longdat$Alcohol > 3, 1 , NA)))
  
}

#################################################

# save vars still with NA (needed only in script 3)
mis.var=names(longdat)[colSums(is.na(longdat[(nrow(Data)+1):nrow(longdat),]))!=0]

mod_imp <- as.mids(longdat)


# # ==========================================================
# # Save 
 saveRDS(mod_imp, file = 'adult.moderator_data_fu_imp.rds')

