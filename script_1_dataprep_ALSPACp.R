# CREDIT:
# some of this code was written by
# Camille Souama: https://github.com/camillesouama/earlycause-tools/tree/main/Amsterdam%20UMC/Meta-analysis%20on%20childhood%20maltreatment%20and%20(comorbid)%20depression%20and%20cardiometabolic%20disease
# Serena Defina: https://github.com/SereDef
# and adjusted, merged, extended by Esther Walton

#################### 0. General settings ####################

# Set working directory, adjust with your own pathway
setwd("...") 

# Upload packages
library('foreign')
library('table1')
library('tidyverse')
library('dplyr')
library('gt')
library('plyr')
library('nnet')
library('pscl')
library('broom')
library('ggplot2')
library('car')
library('mice')
library('gtsummary') 

#################### 1. Read data ####################

Data.orig <- read.spss("...",            
                  use.value.labels=FALSE,
                  to.data.frame=TRUE,
                  max.value.labels=Inf)

#### ALSPAC: check for duplicates
table(duplicated(Data.orig[,c("cidB2957")]))

Data.orig=Data.orig %>%
  distinct(cidB2957, .keep_all = TRUE)

table(duplicated(Data.orig[,c("cidB2957")]))

Data.orig <- subset(Data.orig, mult == 0 | is.na(mult))
Data <- subset(Data.orig[,c("cidB2957","qlet", "mult")], mult == 0 | is.na(mult))


#################### Code all variables #################### 

# 4.1. Code childhood maltreatment variables : 1=yes, 0=no

Data$CM_PA <- as.factor(ifelse(Data.orig$pb466a == 1, 1, ifelse(Data.orig$pb466a == 2, 0, NA)))

Data$CM_EA <- as.factor(ifelse(Data.orig$pb470a == 1, 1, ifelse(Data.orig$pb470a == 2, 0, NA)))

Data$CM_SA <- as.factor(ifelse(Data.orig$pb472a == 1, 1, ifelse(Data.orig$pb472a == 2, 0, NA)))

# check
table(Data.orig$pb466a)
table(Data$CM_PA)

# check
table(Data.orig$pb470a)
table(Data$CM_EA)

# check
table(Data.orig$pb472a)
table(Data$CM_SA)

Data$CM <- as.factor(ifelse(Data$CM_PA == 1 | Data$CM_EA == 1 | Data$CM_SA == 1, 1,
             ifelse(Data$CM_PA == 0 & Data$CM_EA == 0 & Data$CM_SA == 0, 0, NA)))

# check
table(Data$CM)

Data$CM_sev <- as.factor(rowSums(data.frame(as.numeric(as.character(Data$CM_PA)),
                                  as.numeric(as.character(Data$CM_EA)),
                                  as.numeric(as.character(Data$CM_SA)))))

# check
table(Data$CM_sev)

# 2.2. Code depression-related variables : 1=yes, 0=no 

# reminder: fu1=original coding: dx at follow-up, regardless of past tp's
# fu = dx at follow-up, not present at baseline (when LS factors were measured)
# baseline = dx at baseline

Data$MDD.fu1 <- ifelse(Data.orig$fa3255 >= 13, 1, ifelse(Data.orig$fa3255 < 13 & Data.orig$fa3255 >= 0, 0, NA))
Data$MDD.baseline1 <- Data.orig$pb260
Data$MDD.baseline2 <- ifelse(Data.orig$pc102 >= 13, 1, ifelse(Data.orig$pc102 < 13 & Data.orig$pc102 >= 0, 0, NA))
Data$MDD.baseline3 <- ifelse(Data.orig$pd200 >= 13, 1, ifelse(Data.orig$pd200 < 13 & Data.orig$pd200 >= 0, 0, NA))
Data$MDD.baseline4 <- ifelse(Data.orig$pe290 >= 13, 1, ifelse(Data.orig$pe290 < 13 & Data.orig$pe290 >= 0, 0, NA))

var_MDD=c('MDD.baseline1', 'MDD.baseline2', 'MDD.baseline3', 'MDD.baseline4')

# MDD at any baseline
Data$MDD.baseline <- ifelse(rowSums(Data[,var_MDD], na.rm = T)>0 , 1, 
                      ifelse(rowSums(is.na(Data[,var_MDD]))==length(var_MDD), NA, 0))

# MDD at follow-up
Data$MDD.fu <- ifelse(Data$MDD.baseline==0 & Data$MDD.fu1==1, 1, 
                      ifelse(rowSums(Data[,c(var_MDD,'MDD.fu1')], na.rm = T)==0, 0, NA))




# 2.3. Code cardiometabolic-related variables : 1=yes, 0=no

Data$CD <- as.factor(ifelse(Data.orig$fa5000 == 1, 1,
                  ifelse(Data.orig$fa5000 == 2, 0, NA)))

# check
table(Data.orig$fa5000)
table(Data$CD)

Data$CVD1 <- as.factor(ifelse(Data.orig$fa5100 == 1, 1,
                    ifelse(Data.orig$fa5100 == 2, 0, NA)))

# check
table(Data.orig$fa5100)
table(Data$CVD1)

#diabtes scores
#ever had diabetes

# fu1
Data$DM.fu1 <- as.factor(ifelse(Data.orig$fa5230 == 1, 1,
                  ifelse(Data.orig$fa5230 == 2, 0, NA)))

# baseline
Data$DM.baseline <- as.factor(ifelse(Data.orig$pl1070 == 1, 1,
                               ifelse(Data.orig$pl1070 == 2, 0, NA)))

# only at follow-up
Data$DM.fu <- ifelse(Data$DM.baseline==0 & Data$DM.fu1==1, 1, 
                      ifelse(Data$DM.fu1==0, 0, NA))

Data$DM.fu=as.factor(Data$DM.fu)


# check
View(Data[,grep('^DM', names(Data))])
table(Data$DM.fu)

# 4.4. Code PCM multimorbidity variables : 0=healthy, 1=MDD, 2=CMD, 3=comorbid
Data$CMD <- as.factor(ifelse(Data$CVD1 == 0 & Data$DM.fu1 == 0,0,
                   ifelse(Data$CVD1 == 1 | Data$DM.fu1 == 1,1,NA)))

Data$CMD.fu <- as.factor(ifelse(Data$CVD1 == 0 & Data$DM.fu == 0,0,
                             ifelse(Data$CVD1 == 1 | Data$DM.fu == 1,1,NA)))

Data$PCM <- as.factor(ifelse(Data$MDD.fu1 == 0 & Data$CMD == 0, 0,
                   ifelse(Data$MDD.fu1 == 1 & Data$CMD == 0, 1,
                          ifelse(Data$MDD.fu1 == 0 & Data$CMD == 1, 2,
                                 ifelse(Data$MDD.fu1 == 1 & Data$CMD == 1, 3, NA)))))


Data$PCM.fu <- as.factor(ifelse(Data$MDD.fu == 0 & Data$CMD == 0, 0,
                             ifelse(Data$MDD.fu == 1 & Data$CMD == 0, 1,
                                    ifelse(Data$MDD.fu == 0 & Data$CMD == 1, 2,
                                           ifelse(Data$MDD.fu == 1 & Data$CMD == 1, 3, NA)))))

# check
table(Data$PCM)
table(Data$PCM.fu)


# 2.5. Code covariate variables : change - age at followup 
Data$Age <- Data.orig$fa9992

# check
summary(Data.orig$fa9992)
summary(Data$Age)

# only males
Data$Sex <- 0

# recode as factor and check
Data$Sex = as.factor(Data$Sex)
table(Data$Sex)


Data$Educ <- as.factor(ifelse(Data.orig$pb321 == 1, 2,
                    ifelse(Data.orig$pb321 == 2 & Data.orig$pb312 == 1, 1, 
                              ifelse(Data.orig$pb321 == 2 & Data.orig$pb312 == 2 & Data.orig$pb311 == 1 | Data.orig$pb311 == 2 , 0, NA))))

# check
table(Data.orig$pb311) # O-level 1=y
table(Data.orig$pb312) # A-level 1=y
table(Data.orig$pb321) # uni 1=y
table(Data$Educ) #  ok

# ethnicity code adapted from: https://github.com/SereDef/association-ELS-PCM-project/blob/main/ALSPAC/3.1.PCM_outcomes_covs_aux.R
# should be ethn; 0 = European ancestry, 1 = non-European ancestry
# in ALSPAC: 1 = white, >1 = non-white
Data$Ethn <- ifelse(is.na(Data.orig$c801), NA, 
                    ifelse(Data.orig$c801 == 1 , 0, 1))

# recode as factor and check
Data$Ethn = as.factor(Data$Ethn)
table(Data.orig$c801)
table(Data$Ethn)


# exercise at 7y follow-up
# 1	Never
# 2	<= once a month
# 3	<= once a week
# 4	2-3 times a week
# 5	4-5 times a week
# 6	Most days

# table(Data.orig$pk5230) #Frequency during the past year partner spent time hiking or walking
# table(Data.orig$pk5231) #Frequency during the past year partner spent time jogging
# table(Data.orig$pk5232) #Frequency during the past year partner spent time running
# table(Data.orig$pk5233) #Frequency during the past year partner spent time cycling
# table(Data.orig$pk5234) #Frequency during the past year partner spent time doing keep fit, aerobics, step aerobics, etc
# table(Data.orig$pk5235) #Frequency during the past year partner spent time playing tennis, squash, badminton, etc
# table(Data.orig$pk5236) #Frequency during the past year partner spent time swimming
# table(Data.orig$pk5237) #Frequency during the past year partner spent time doing another energetic leisure activity

var_phyact=c('pk5230','pk5231','pk5232', 'pk5233', 'pk5234', 'pk5235', 'pk5236', 'pk5237')

Data$PhyAct=apply(Data.orig[,var_phyact],1, max, na.rm=T)
Data$PhyAct[is.infinite(Data$PhyAct)] <- NA




# following WHO guidelines of >3 times p/w 

Data$PhyAct_bin = as.factor(ifelse(Data$PhyAct>4,1,
                                   ifelse(Data$PhyAct<5,0,NA)))


table(Data$PhyAct, Data$PhyAct_bin)

#Cigarettes smoked per day (larger=bad)

# #pd620
# 0	None
# 1	1-4	
# 5	5-9	
# 10	10-14	
# 15	15-19
# 20	20-24
# 25	25-29	
# 30	30+
# 
# #pe450
# 0	None	
# 1	1-4	
# 5	5-9	
# 8	Pipe Only	EW => move to 1
# 9	Cigars Only	EW => move to 1
# 10	10-14	
# 15	15-19	
# 20	20-24	
# 25	25-29	
# 30	30+
# 
Data.orig$pe450 = car::recode(Data.orig$pe450, "8 = 1; 9 = 1") 
# 
# #pf7090
# 0	None
# 1	1-4	
# 5	5-9	
# 8	Pipe only EW => move to 1
# 9	Cigars only	EW => move to 1
# 10	10-14	
# 15	15-19	
# 20	20-24	
# 25	25-29	
# 30	30+	
# 
Data.orig$pf7090 = car::recode(Data.orig$pf7090, "8 = 1; 9 = 1") 
# 
# #ph6180
# 0	None	
# 1	1-4	
# 5	5-9	
# 8	Pipe only	EW => move to 1
# 9	Cigars only	EW => move to 1
# 10	10-14	
# 15	15-19	
# 20	20-24	
# 25	25-29	
# 30	30+	
# 97	Occasionally EW => move to 1	
Data.orig$ph6180 = car::recode(Data.orig$ph6180, "8 = 1; 9 = 1; 97=1") 


var_smok=c('pd620','pe450','pf7090','ph6180')

Data$Smoke <- rowMeans(Data.orig[,var_smok],na.rm=T)


Data$Smoke_bin <- as.factor(ifelse(Data$Smoke>0,1,
                               ifelse(Data$Smoke==0,0,NA)))

table(Data$Smoke,Data$Smoke_bin)

# alcohol (larger=bad)
# 1	Never drink alcohol	
# 2	Very occasionally (less than once a week)	
# 3	Occasionally (at least once a week)	
# 4	Drink 1-2 glasses nearly every day	
# 5	Drink 3-9 glasses every day	
# 6	Drink at least 10 glasses a day

Mode <- function(x, na.rm = TRUE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

var_alc=c('pe452', 'pf7100', 'pg5050', 'ph6190')

Data$Alcohol <- apply(Data.orig[,var_alc],1,Mode)


# binary cutoff: NHS
# it's recommended to drink no more than 14 units of alcohol a week, spread across 3 days or more. That's around 6 medium (175ml) glasses of wine, or 6 pints of 4% beer.
# so, cat=4 ok, but 5 or 6 isn't
Data$Alcohol_bin <- as.factor(ifelse(Data$Alcohol>4,1,
                                   ifelse(Data$Alcohol<5,0,NA)))

table(Data$Alcohol,Data$Alcohol_bin)



# diet
Data$Diet_healthy=Data.orig$pg2540
Data$Diet_trad=Data.orig$pg2541
Data$Diet_proc.conf=Data.orig$pg2542
Data$Diet_semiveget=Data.orig$pg2543


# median split
Data$Diet_healthy_bin=as.factor(ifelse(Data$Diet_healthy>median(Data$Diet_healthy, na.rm=T),1,
                             ifelse(Data$Diet_healthy<=median(Data$Diet_healthy, na.rm=T),0,NA)))

Data$Diet_trad_bin=as.factor(ifelse(Data$Diet_trad>median(Data$Diet_trad, na.rm=T),1,
                                    ifelse(Data$Diet_trad<=median(Data$Diet_trad, na.rm=T),0,NA)))

Data$Diet_proc.conf_bin=as.factor(ifelse(Data$Diet_proc.conf>median(Data$Diet_proc.conf, na.rm=T),1,
                                         ifelse(Data$Diet_proc.conf<=median(Data$Diet_proc.conf, na.rm=T),0,NA)))

Data$Diet_semiveget_bin=as.factor(ifelse(Data$Diet_semiveget>median(Data$Diet_semiveget, na.rm=T),1,
                                         ifelse(Data$Diet_semiveget<=median(Data$Diet_semiveget, na.rm=T),0,NA)))


## remove those with less than 50% data
percent_missing <- function(var) { sum(is.na(var)) / length(var) * 100 }
Data$percent_missing <- apply(Data[,2:ncol(Data)], 1, percent_missing)

## add descriptives table here
col_order <- c("Age", "Sex", "Ethn",'Educ',
               "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
               'MDD.fu1','DM.fu1','CD','CMD','CVD1', 'PCM',
               'Alcohol','Alcohol_bin',
               'Smoke','Smoke_bin',
               'PhyAct','PhyAct_bin')

## add follow-up descriptives table here
col_order.fu <- c("Age", "Sex", "Ethn",'Educ',
               "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
               'MDD.fu','DM.fu1','CD','CMD','CVD1', 'PCM.fu',
               'Alcohol','Alcohol_bin',
               'Smoke','Smoke_bin',
               'PhyAct','PhyAct_bin')

Data1 <- Data[, c('cidB2957',col_order, grep("Diet_h", names(Data), value=T),
                  grep("Diet_t", names(Data), value=T),
                 grep("veg", names(Data), value=T),
                 grep("Diet_con", names(Data), value=T),
                 grep("Diet_pro", names(Data), value=T),'percent_missing')]

Data3 <- Data[, c('cidB2957',col_order.fu, grep("Diet_h", names(Data), value=T),
                  grep("Diet_t", names(Data), value=T),
                  grep("veg", names(Data), value=T),
                  grep("Diet_con", names(Data), value=T),
                  grep("Diet_pro", names(Data), value=T),'percent_missing')]


Data1 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)", missing = 'always') %>%
  modify_caption("**Original Sample (full, not imputed) Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = "...")#, pr_section = sect_properties)

data.cor <- data.frame(lapply(Data1 , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="...")


Data2 = Data1[Data1$percent_missing<50,]

Data2 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)", missing = 'always') %>%
  modify_caption("**Original Sample (50% nonmiss, not imputed) Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = "...")#, pr_section = sect_properties)

data.cor <- data.frame(lapply(Data2 , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="...")

rm(list=c('Data.orig', 'c', 'data.cor'))

# follow-up
Data3 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)",missing = 'always') %>%
  modify_caption("**Original Sample at follow-up (full, not imputed) Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = "...")#, pr_section = sect_properties)

data.cor <- data.frame(lapply(Data3 , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="./results/Sample_orig.fu.full.corr.csv")


Data4 = Data3[Data3$percent_missing<50,]

Data4 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)",missing = 'always') %>%
  modify_caption("**Original Sample at follow-up (50% nonmiss, not imputed) Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = "...")#, pr_section = sect_properties)

data.cor <- data.frame(lapply(Data4 , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="...")

Data = Data[Data$percent_missing<50,]

# optional: save Main file
# saveRDS(Data, file = '...')


