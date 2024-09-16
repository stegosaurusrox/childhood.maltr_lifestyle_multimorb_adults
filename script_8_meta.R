library('ggmice')
library('ggplot2')
library('ggsci')
library('forestplot')
library('patchwork')
library('openxlsx')
library('metafor')

setwd("...")


# define, which mod (bin or cont) you want to plot

mod_sheet='CM_modz_PCM_int' # CM_modz_PCM_int or CM_modb_PCM_int (earlier versions also without interaction: CM_modz_PCM_noint, CM_modb_PCM_noint)
mod_type='_z' # _z or _b

########### partner results

results_plot_PCM_ALSPACp <- read.xlsx('...', sheet = mod_sheet)
results_plot_PCM_ALSPACp <- results_plot_PCM_ALSPACp[,c('model', 'y.level', 'term', 'estimate', 'std.error', 'p.value', 'sign', 'OR', 'lci', 'uci')]
results_plot_PCM_ALSPACp$model <- paste0(results_plot_PCM_ALSPACp$model,'_ALSPAC.partners')
results_plot_PCM_ALSPACp$dataset <- 'ALSPAC partners'


############ UKBB results

results_plot_PCM_UKBB <- read.xlsx('...', sheet = mod_sheet)
results_plot_PCM_UKBB <- results_plot_PCM_UKBB[,c('model', 'y.level', 'term', 'estimate', 'std.error', 'p.value', 'sign', 'OR', 'lci', 'uci')]
results_plot_PCM_UKBB$model <- paste0(results_plot_PCM_UKBB$model,'_UKBB')
results_plot_PCM_UKBB$dataset <- 'UK Biobank'

########## ALSPAC mother results

results_plot_PCM_ALSPACm <- read.xlsx('...', sheet = mod_sheet)
results_plot_PCM_ALSPACm <- results_plot_PCM_ALSPACm[,c('model', 'y.level', 'term', 'estimate', 'std.error', 'p.value', 'sign', 'OR', 'lci', 'uci')]
results_plot_PCM_ALSPACm$model <- paste0(results_plot_PCM_ALSPACm$model,'_ALSPAC.mothers')
results_plot_PCM_ALSPACm$dataset <- 'ALSPAC mothers'

########## NESDA results

results_plot_PCM_NESDA <- read.xlsx('...', sheet = mod_sheet)
results_plot_PCM_NESDA <- results_plot_PCM_NESDA[,c('model', 'y.level', 'term', 'estimate', 'std.error', 'p.value', 'sign', 'OR', 'lci', 'uci')]
results_plot_PCM_NESDA$model <- paste0(results_plot_PCM_NESDA$model,'_NESDA')
results_plot_PCM_NESDA$dataset <- 'NESDA'

results_plot_PCM=rbind(results_plot_PCM_ALSPACm,results_plot_PCM_ALSPACp, results_plot_PCM_NESDA, results_plot_PCM_UKBB)
results_plot_PCM=results_plot_PCM[order(results_plot_PCM$model),]

##### extract, rename, and forest plot
# 
results_plot_PCM$y.level <- gsub("intern", "MDD", results_plot_PCM$y.level)
results_plot_PCM$y.level <- gsub("adipos", "CMD", results_plot_PCM$y.level)

head(results_plot_PCM)

## optional: remove diet
results_plot_PCM=results_plot_PCM[grep('Diet', results_plot_PCM$model, invert = T),]

## optional: only diet and order
# results_plot_PCM=results_plot_PCM[c(grep('Diet_h', results_plot_PCM$model),
#                                     grep('Diet_v|Diet_semi', results_plot_PCM$model),
#                                     grep('Diet_conf|Diet_proc_|Diet_proc.conf', results_plot_PCM$model)),]

########################### metafor

# export as pdf A4 landscape
par(mfrow = c(3, 3),mar = c(2,2,2,2))

# CM on MDD (alcohol model)
Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1$', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1



sum_forest_plots=function(data=data, title){
  RE_model1 <- rma(yi = estimate , sei = std.error, data = data,
                   slab=dataset)
  summary(RE_model1)
  
  OR_model1 <- predict(RE_model1, transf=exp, digits=3)
  weights_model1 <- paste0(formatC(weights(RE_model1), format="f", digits=1, width=4), "%")
  
  forest(RE_model1, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,10),
         ilab=cbind(NA, weights_model1), ilab.xpos=c(-9.5,-4.5),
         cex=1.10, 
         header=paste(title,'\n',"Study (weight)"), mlab="RE Model", xlab=' ',
         shade=TRUE)
  
  ## add text with Q-value, dfs, p-value, and I^2 statistic
  text(-10, -1, pos=4, #cex=0.75,
       bquote(paste("(",I^2, " = ", .(fmtx(RE_model1$I2, digits=1)), "%)")))
  
}


sum_forest_plots(data=Model1, title="CM on Depression (alcohol model)")

##############


# CM on MDD (physical activity model)
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1$', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Depression (physical activity model)")




##############


# CM on MDD (smoking model)
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('CM1$', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Depression (smoking model)")



################ 

# CM on CMD (alcohol model)

Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1$', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Cardiometabolic Disease (alcohol model)")


##############


# CM on CMD (physical activity model)
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1$', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Cardiometabolic Disease (physical activity model)")


##############


# CM on CMD (smoking model)
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('CM1$', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Cardiometabolic Disease (smoking model)")


################ 

#CM on PCM (alcohol model)

Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1$', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Comorbidity (alcohol model)")


##############


# CM on PCM (physical activity model)
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1$', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Comorbidity (physical activity model)")


##############


# CM on PCM (smoking model)
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('CM1$', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM on Comorbidity (smoking model)")


######################## alcohol on MDD

# export as pdf A4 landscape
par(mfrow = c(3, 3),mar = c(2,2,2,2))

Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^Alcohol_z$', Model1$term) ,]
#Model1=Model1[grep('^Alcohol_bin1$', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Alcohol on Depression")


##############


# Physical Activity on MDD

Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^PhyAct_z$', Model1$term) ,]
#Model1=Model1[grep('^PhyAct_bin*', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Physical Activity on Depression")


##############


# smoking on MDD
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('^Smoke_z$', Model1$term) ,]
#Model1=Model1[grep('^Smoke_bin*', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Smoking on Depression")


################ 

# alcohol on CMD
Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^Alcohol_z$', Model1$term) ,]
#Model1=Model1[grep('^Alcohol_bin*', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Alcohol on Cardiometabolic Disease")


##############


# Physical Activity on CMD
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^PhyAct_z$', Model1$term) ,]
#Model1=Model1[grep('^PhyAct_bin*', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Physical Activity on Cardiometabolic Disease")

##############


# smoking on CMD
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('^Smoke_z$', Model1$term) ,]
#Model1=Model1[grep('^Smoke_bin*', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Smoking on Cardiometabolic Disease")

################ 

# alcohol on PCM
Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^Alcohol_z$', Model1$term) ,]
#Model1=Model1[grep('^Alcohol_bin*', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Alcohol on Comorbidity")


##############


# Physical Activity on PCM 
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^PhyAct_z$', Model1$term) ,]
#Model1=Model1[grep('^PhyAct_bin*', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Physical Activity on Comorbidity")

##############


# smoking on PCM
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('^Smoke_z$', Model1$term) ,]
#Model1=Model1[grep('^Smoke_bin*', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="Smoking on Comorbidity")

######################## interactions

######################## CM*alcohol on MDD

# export as pdf A4 landscape
par(mfrow = c(3, 3),mar = c(2,2,2,2))

Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:Alcohol_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:Alcohol_bin*', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Alcohol on Depression")


##############


# CM*Physical Activity on MDD

Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:PhyAct_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:PhyAct_bin*', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Physical Activity on Depression")


##############


# CM*smoking on MDD
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:Smoke_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:Smoke_bin*', Model1$term) ,]
Model1=Model1[grep('MDD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Smoking on Depression")


################ 

# CM*alcohol on CMD
Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:Alcohol_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:Alcohol_bin*', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Alcohol on Cardiometabolic Disease")


##############


# CM*Physical Activity on CMD
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:PhyAct_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:PhyAct_bin*', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Physical Activity on Cardiometabolic Disease")

##############


# CM*smoking on CMD
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:Smoke_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:Smoke_bin*', Model1$term) ,]
Model1=Model1[grep('CMD', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Smoking on Cardiometabolic Disease")

################ 

# CM*alcohol on PCM
Model1=results_plot_PCM[grep('^Alcohol', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:Alcohol_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:Alcohol_bin*', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Alcohol on Comorbidity")


##############


# CM*Physical Activity on PCM 
Model1=results_plot_PCM[grep('^PhyAct', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:PhyAct_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:PhyAct_bin*', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Physical Activity on Comorbidity")

##############


# CM*smoking on PCM
Model1=results_plot_PCM[grep('^Smoke', results_plot_PCM$model) ,]
Model1=Model1[grep('^CM1:Smoke_z$', Model1$term) ,]
#Model1=Model1[grep('^CM1:Smoke_bin*', Model1$term) ,]
Model1=Model1[grep('comorb', Model1$y.level) ,]
Model1

sum_forest_plots(data=Model1, title="CM * Smoking on Comorbidity")


