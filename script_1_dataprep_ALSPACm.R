# CREDIT:
# some of this code was written by
# Camille Souama: https://github.com/camillesouama/earlycause-tools/tree/main/Amsterdam%20UMC/Meta-analysis%20on%20childhood%20maltreatment%20and%20(comorbid)%20depression%20and%20cardiometabolic%20disease
# Serena Defina: https://github.com/SereDef
# and adjusted, merged, extended by Esther Walton

#################### 0. General settings ####################

# Set working directory, adjust with your own pathway
setwd("...") 

#install packages below
install.packages(foreign)


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

# 2.1. Code childhood maltreatment variables

Data$CM_PA <- as.factor(ifelse(Data.orig$c416a == 1, 1, ifelse(Data.orig$c416a == 2, 0, NA)))

Data$CM_EA <- as.factor(ifelse(Data.orig$c420a == 1, 1, ifelse(Data.orig$c420a == 2, 0, NA)))

Data$CM_SA <- as.factor(ifelse(Data.orig$c422a == 1, 1, ifelse(Data.orig$c422a == 2, 0, NA)))

# check
table(Data.orig$c416a)
table(Data$CM_PA)

# check
table(Data.orig$c420a)
table(Data$CM_EA)

# check
table(Data.orig$c422a)
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

# 2.2. Code depression-related variables: 1=yes, 0=no 
 
# reminder: fu1=original coding: dx at follow-up, regardless of past tp's
# fu = dx at follow-up, not present at baseline (when LS factors were measured)
# baseline = dx at baseline

Data$MDD.fu1 <- Data.orig$t3255

# careful: b370 and c600 are coded as 0/1 despite the documentation saying it's the EDPS score (0-30)
Data$MDD.baseline1 <- ifelse(Data.orig$b370 == 1, 1, ifelse(Data.orig$b370 == 0, 0, NA))
Data$MDD.baseline2 <- ifelse(Data.orig$c600 == 1, 1, ifelse(Data.orig$c600 == 0, 0, NA))
Data$MDD.baseline3 <- ifelse(Data.orig$e390 >= 13, 1, ifelse(Data.orig$e390 < 13 & Data.orig$e390 >= 0, 0, NA))
Data$MDD.baseline4 <- ifelse(Data.orig$f200 >= 13, 1, ifelse(Data.orig$f200 < 13 & Data.orig$f200 >= 0, 0, NA))
Data$MDD.baseline5 <- ifelse(Data.orig$g290  >= 13, 1, ifelse(Data.orig$g290  < 13 & Data.orig$g290  >= 0, 0, NA))



var_MDD=c('MDD.baseline1', 'MDD.baseline2', 'MDD.baseline3', 'MDD.baseline4', 'MDD.baseline5')

# MDD at any baseline
Data$MDD.baseline <- ifelse(rowSums(Data[,var_MDD], na.rm = T)>0 , 1, 
                            ifelse(rowSums(is.na(Data[,var_MDD]))==length(var_MDD), NA, 0))

# MDD at follow-up
Data$MDD.fu <- ifelse(Data$MDD.baseline==0 & Data$MDD.fu1==1, 1,
                      ifelse(rowSums(Data[,c(var_MDD,'MDD.fu1')], na.rm = T)==0, 0, NA))

# 2.3. Code cardiometabolic-related variables : 1=yes, 0=no

Data$CD <- as.factor(ifelse(Data.orig$t5000 == 1, 1,
                            ifelse(Data.orig$t5000 == 2, 0, NA)))

# check
table(Data.orig$t5000)
table(Data$CD)

Data$CVD1 <- as.factor(ifelse(Data$CD == 1 | Data.orig$t5100 == 1, 1,
                              ifelse(Data$CD == 0 & Data.orig$t5100 == 2, 0, NA)))

# check
table(Data.orig$t5100)
table(Data$CVD1) # diff, but that's ok, as also included CD

Data$CVD2 <- as.factor(ifelse(Data$CVD1 == 1 | Data.orig$t5060 == 1, 1,
                              ifelse(Data$CVD1 == 0 & 
                                       Data.orig$t5060 == 2, 0, NA)))

table(Data.orig$t5060)
table(Data$CVD2) # diff, but that's ok, as also included CVD1


#medication variables not used in ALSPAC

#diabtes scores
#ever had diabetes
Data$DM <- as.factor(ifelse(Data.orig$t5230 == 1, 1,
                            ifelse(Data.orig$t5230 == 2, 0, NA)))

# check
table(Data.orig$t5230)
table(Data$DM)


# 2.4. Code PCM multimorbidity variables 
Data$CMD <- as.factor(ifelse(Data$CVD1 == 0 & Data$DM == 0,0,
                             ifelse(Data$CVD1 == 1 | Data$DM == 1,1,NA)))

Data$CardiacMD <- as.factor(ifelse(Data$CD == 0 & Data$DM == 0,0,
                                   ifelse(Data$CD == 1 | Data$DM == 1,1,NA)))

Data$CVD2MD <- as.factor(ifelse(Data$CVD2 == 0 & Data$DM == 0,0,
                                ifelse(Data$CVD2 == 1 | Data$DM == 1,1,NA)))


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


# 2.5. Code covariate variables: change - age at followup 
Data$Age <- Data.orig$t9994

# check
summary(Data.orig$t9994)
summary(Data$Age)

# and age at baseline (only for descritives)
Data$Age.baseline <- Data.orig$mz028b
summary(Data.orig$mz028b)
summary(Data$Age.baseline)

# only females
Data$Sex <- 1

# recode as factor and check
Data$Sex = as.factor(Data$Sex)
table(Data$Sex)


Data$Educ <- as.factor(ifelse(Data.orig$c645a == 1 | Data.orig$c645a == 2 | Data.orig$c645a == 3 | Data.orig$c645a == 0, 0,
                              ifelse(Data.orig$c645a == 4, 1, ifelse(Data.orig$c645a == 5, 2, NA))))
# check
table(Data.orig$c645a) 
table(Data$Educ) 

# ethnicity code adapted from: https://github.com/SereDef/association-ELS-PCM-project/blob/main/ALSPAC/3.1.PCM_outcomes_covs_aux.R
# should be ethn; 0 = European ancestry, 1 = non-European ancestry
# in ALSPAC: 1 = white, >1 = non-white
Data$Ethn <- ifelse(is.na(Data.orig$c800), NA, 
                    ifelse(Data.orig$c800 == 1 , 0, 1))

# recode as factor and check
Data$Ethn = as.factor(Data$Ethn)
table(Data.orig$c800)
table(Data$Ethn)

############################# exercise at 7y, 8y pr 11y follow-up
# more frequ per tp, then standardizem and then across tp's: mean and rank


############ 7y

table(Data.orig$m5236) #Average time spent swimming in past year
# 1	Never
# 2	Once a month or less
# 3	Once a week or less
# 4	2-3 times a week
# 5	4-5 times a week
# 6	Most days

table(Data.orig$m5232) #Average time spent running in past year
# 1	Never
# 2	Once a month or less
# 3	Once a week or less
# 4	2-3 times a week
# 5	4-5 times a week
# 6	Most days

table(Data.orig$m5233) #Average time spent cycling in past year
# 1	Never
# 2	Once a month or less
# 3	Once a week or less
# 4	2-3 times a week
# 5	4-5 times a week
# 6	Most days

var_phyact.m=c('m5236','m5232','m5233')

Data$PhyAct.m=apply(Data.orig[,var_phyact.m],1, max, na.rm=T)
Data$PhyAct.m[is.infinite(Data$PhyAct.m)] <- NA

#View(Data.orig[,c(var_phyact.m,'PhyAct.m')])
table(Data$PhyAct.m) # 5 or 6 = good

Data$PhyAct.m.std=scale(Data$PhyAct.m)
Data$PhyAct.m.bin=ifelse(Data$PhyAct.m>4,1,
                                   ifelse(Data$PhyAct.m<5,0,NA))


table(Data$PhyAct.m.bin)

################# 8y

table(Data.orig$n5119) #Amount of time mother spends cycling per week
# 1	>6 hours
# 2	2 - 6 hours
# 3	<2 hours
# 4	None

table(Data.orig$n5115) #Amount of time mother spends playing tennis or badminton per week
# 1	>6 hours
# 2	2 - 6 hours
# 3	<2 hours
# 4	None

table(Data.orig$n5111) #Amount of time mother spends doing aerobics per week
# 1	>6 hours
# 2	2 - 6 hours
# 3	<2 hours
# 4	None

table(Data.orig$n5116) #Amount of time mother spends swimming per week
# 1	>6 hours
# 2	2 - 6 hours
# 3	<2 hours
# 4	None

var_phyact.n=c('n5119','n5115','n5111', 'n5116')

Data$PhyAct.n=apply(Data.orig[,var_phyact.n],1, min, na.rm=T)
Data$PhyAct.n[is.infinite(Data$PhyAct.n)] <- NA

#View(Data[,c(var_phyact.n,'PhyAct.n')])
table(Data$PhyAct.n) # 1 or 2 = good

# reverse code
Data$PhyAct.n.rev = Data$PhyAct.n*(-1) # 1 or 2 = good
table(Data$PhyAct.n.rev) # -1 or -2 = good
Data$PhyAct.n.std=scale(Data$PhyAct.n.rev)

# WHO recommends at least 150 minutes of moderate physical activity, or 75 minutes of vigorous physical activity per week => recode to 1=good if at least 2h, 0=bad if less
# so, 2 (=2-6h) or 1 = good

Data$PhyAct.n.bin=ifelse(Data$PhyAct.n<3,1,
                         ifelse(Data$PhyAct.n>2,0,NA))

table(Data$PhyAct.n.bin)


################# 11y

table(Data.orig$r8201) #Frequency respondent played tennis/badminton in past year
Data.orig$r8201 = car::recode(Data.orig$r8201, "0 = NA") 
# 1	Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

table(Data.orig$r8271) #Frequency respondent played netball, volleyball, basketball in past year
Data.orig$r8271 = car::recode(Data.orig$r8271, "0 = NA") 
# 1 Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

table(Data.orig$r8241) #Frequency respondent played football/hockey in past year
Data.orig$r8241 = car::recode(Data.orig$r8241, "0 = NA") 
# 1 Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

table(Data.orig$r8051) #Frequency respondent has gone cycling for pleasure in past year
Data.orig$r8051 = car::recode(Data.orig$r8051, "0 = NA") 
# 1 Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

table(Data.orig$r8211) #Frequency respondent played squash in past year
Data.orig$r8211 = car::recode(Data.orig$r8211, "0 = NA") 
# 1 Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

table(Data.orig$r8111) #Frequency respondent did high impact, step aerobics in past year
Data.orig$r8111 = car::recode(Data.orig$r8111, "0 = NA") 
# 1 Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

table(Data.orig$r8121) #Frequency respondent did other types of aerobics in past year
Data.orig$r8121 = car::recode(Data.orig$r8121, "0 = NA") 
# 1 Every day
# 2	3-6 times a week
# 3	Once or twice a week
# 4	1-3 times a month
# 5	Less than once a month
# 6	None

var_phyact.r=c('r8201','r8241','r8211', 'r8051', 'r8111', 'r8121')

Data$PhyAct.r=apply(Data.orig[,var_phyact.r],1, min, na.rm=T)
Data$PhyAct.r[is.infinite(Data$PhyAct.r)] <- NA

table(Data$PhyAct.r) # 1 or 2 = good

# reverse code
Data$PhyAct.r.rev=Data$PhyAct.r*(-1)
table(Data$PhyAct.r.rev) # -1 or -2 = good

Data$PhyAct.r.std=scale(Data$PhyAct.r.rev)
Data$PhyAct.r.bin=ifelse(Data$PhyAct.r<3,1,
                         ifelse(Data$PhyAct.r>2,0,NA))


table(Data$PhyAct.r.bin)

# mean across std m/n/r values
Data$PhyAct.std=apply(Data[,c('PhyAct.n.std','PhyAct.m.std', 'PhyAct.r.std')],1, mean, na.rm=T)
Data$PhyAct.std[Data$PhyAct.std == 'NaN'] <- NA
Data$PhyAct.rank=as.numeric(as.factor(Data$PhyAct.std))
Data$PhyAct=Data$PhyAct.rank


# following WHO guidelines of >3 times p/w at at least 1 of the 3 tps
Data$PhyAct_bin = as.factor(ifelse(rowSums(Data[,c('PhyAct.m.bin', 'PhyAct.n.bin', 'PhyAct.r.bin')], na.rm = T)>0,1,
                                   ifelse(rowSums(Data[,c('PhyAct.m.bin', 'PhyAct.n.bin', 'PhyAct.r.bin')], na.rm = T)<1 &
                                            rowSums(is.na(Data[,c('PhyAct.m.bin', 'PhyAct.n.bin', 'PhyAct.r.bin')]))!=length(c('PhyAct.m.bin', 'PhyAct.n.bin', 'PhyAct.r.bin')),0,NA)))


boxplot(Data$PhyAct_bin, Data$PhyAct)
table(Data$PhyAct_bin) 

# clean up
Data=Data[ , -which(names(Data) %in% c("PhyAct.m","PhyAct.n",'PhyAct.r',
                                  'PhyAct.m.std', 'PhyAct.n.std', 'PhyAct.r.std',
                                  'PhyAct.m.bin', 'PhyAct.n.bin', 'PhyAct.r.bin',
                                  'PhyAct.n.rev', 'PhyAct.r.rev',
                                  'PhyAct.rank', 'PhyAct.std'))]

######################### smoking

#Cigarettes smoked per day (larger=bad)

# recode cigars and occ to 1
# #g820
# #h720
#k6180
Data.orig$k6180 = car::recode(Data.orig$k6180, "9 = 1; 97 = 1") 

#m5160
Data.orig$m5160 = car::recode(Data.orig$m5160, "9 = 1; 97 = 1") 


var_smok=c('g820','h720','k6180','m5160')

Data$Smoke <- rowMeans(Data.orig[,var_smok],na.rm=T)
Data$Smoke[Data$Smoke == 'NaN'] <- NA


Data$Smoke_bin <- as.factor(ifelse(Data$Smoke>0,1,
                                   ifelse(Data$Smoke==0,0,NA)))

table(Data$Smoke,Data$Smoke_bin)


################# alcohol

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

var_alc=c('g822', 'h723', 'k6190')

Data$Alcohol <- apply(Data.orig[,var_alc],1,Mode)


# binary cutoff: NHS
# it's recommended to drink no more than 14 units of alcohol a week, spread across 3 days or more. That's around 6 medium (175ml) glasses of wine, or 6 pints of 4% beer.
# so, cat=4 ok, but 5 or 6 isn't
Data$Alcohol_bin <- as.factor(ifelse(Data$Alcohol>4,1,
                                     ifelse(Data$Alcohol<5,0,NA)))

table(Data$Alcohol,Data$Alcohol_bin)

########################################## diet

Data$Diet_healthy=Data.orig$j1100
Data$Diet_proc=Data.orig$j1101
Data$Diet_conf=Data.orig$j1102
Data$Diet_veget=Data.orig$j1103

## ALSPAC only:
# ref: Northstone, K., & Emmett, P. (2008). A comparison of methods to assess changes in dietary patterns from pregnancy to 4 years post-partum obtained using principal components analysis. British Journal of Nutrition, 99(5), 1099-1106. doi:10.1017/S0007114507842802

# median split
Data$Diet_healthy_bin=as.factor(ifelse(Data$Diet_healthy>median(Data$Diet_healthy, na.rm=T),1,
                             ifelse(Data$Diet_healthy<=median(Data$Diet_healthy, na.rm=T),0,NA)))

Data$Diet_proc_bin=as.factor(ifelse(Data$Diet_proc>median(Data$Diet_proc, na.rm=T),1,
                                    ifelse(Data$Diet_proc<=median(Data$Diet_proc, na.rm=T),0,NA)))

Data$Diet_conf_bin=as.factor(ifelse(Data$Diet_conf>median(Data$Diet_conf, na.rm=T),1,
                                         ifelse(Data$Diet_conf<=median(Data$Diet_conf, na.rm=T),0,NA)))

Data$Diet_veget_bin=as.factor(ifelse(Data$Diet_veget>median(Data$Diet_veget, na.rm=T),1,
                                         ifelse(Data$Diet_veget<=median(Data$Diet_veget, na.rm=T),0,NA)))


## remove those with less than 50% data
percent_missing <- function(var) { sum(is.na(var)) / length(var) * 100 }
Data$percent_missing <- apply(Data[,2:ncol(Data)], 1, percent_missing)


col_order.fu <- c("Age.baseline","Age", "Sex", "Ethn",'Educ',
                  "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
                  'MDD.baseline','MDD.fu1','MDD.fu',
                  'DM','CD','CMD','CVD1', 'PCM','PCM.fu',
                  'Alcohol','Alcohol_bin',
                  'Smoke','Smoke_bin',
                  'PhyAct','PhyAct_bin')


Data3 <- Data[, c('cidB2957',col_order.fu, grep("Diet_h", names(Data), value=T),
                  grep("veg", names(Data), value=T),
                  grep("Diet_con", names(Data), value=T),
                  grep("Diet_pro", names(Data), value=T),'percent_missing')]

Data3$MDD.baseline=as.factor(Data3$MDD.baseline)
Data3$MDD.fu=as.factor(Data3$MDD.fu)
Data3$MDD.fu1=as.factor(Data3$MDD.fu1)



# follow-up
Data3 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)",missing = 'always') %>%
  modify_caption("**Original Sample at baseline & follow-up (full, not imputed) Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = "...")#, pr_section = sect_properties)

data.cor <- data.frame(lapply(Data3 , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="...")


Data4 = Data3[Data3$percent_missing<50,]

Data4 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)",missing = 'always') %>%
  modify_caption("**Original Sample at baseline & follow-up (50% nonmiss, not imputed) Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = "...")#, pr_section = sect_properties)

data.cor <- data.frame(lapply(Data4 , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="./results/Sample_orig.bl.fu.50pc.corr.csv")

Data = Data[Data$percent_missing<50,]

# optional: save Main file
# saveRDS(Data, file = '...')

rm(data.cor,Data.orig, Data3, Data4)
