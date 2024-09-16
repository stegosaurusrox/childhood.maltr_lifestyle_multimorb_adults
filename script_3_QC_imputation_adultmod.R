# CREDIT:
# some of this code was written by
# Serena Defina: https://github.com/SereDef
# and adjusted, merged, extended by Esther Walton

##### prep #####

library(mice)
library(miceadds)


##### load Moderation imputation file 
mod_imp <- readRDS(file = './data/adult.moderator_data_imp.rds')


# Date (for marking output file)
date <- format(Sys.Date(), '%d%m%y')

# Identify the variables which were actually imputed at this round (the moderators) and which don't have NAs (possible as manually set for those you are ill at both time points)
imputed <- names(mod_imp$data[, colSums(is.na(mod_imp$data)) != 0])

# optional: excluse vars with NAs
#imputed <- imputed[imputed!=mis.var]

# Select only imputed variables (this uses miceadds to subset and trnasform back and forth)
mod_imp_qc <- datlist2mids( subset_datlist( mids2datlist(mod_imp), select=c('cidB2957',imputed) ))


# UKBB
mod_imp_qc <- datlist2mids( subset_datlist( mids2datlist(mod_imp), select=c('ID',imputed, mis.var) ))


# Density plot does not work very well with for loops so little trick:

pdf(file = paste0('./results/Density_imp-vs-obs_',date,'.pdf'))


densityplot(mod_imp_qc, ~ CM_PA )
densityplot(mod_imp_qc, ~ CM_EA )
densityplot(mod_imp_qc, ~ CM_SA )
densityplot(mod_imp_qc, ~ CM )
densityplot(mod_imp_qc, ~ CM_sev )
densityplot(mod_imp_qc, ~ MDD.baseline )
densityplot(mod_imp_qc, ~ MDD )
densityplot(mod_imp_qc, ~ AntidepMed)
densityplot(mod_imp_qc, ~ MDD.fu )
densityplot(mod_imp_qc, ~ MDD.fu1 )
densityplot(mod_imp_qc, ~ CVD1 )
densityplot(mod_imp_qc, ~ CVD1.fu )
densityplot(mod_imp_qc, ~ DM )
densityplot(mod_imp_qc, ~ DM_baseline )
densityplot(mod_imp_qc, ~ DM.fu )
densityplot(mod_imp_qc, ~ CMD )
densityplot(mod_imp_qc, ~ PCM )
densityplot(mod_imp_qc, ~ PCM.fu )
densityplot(mod_imp_qc, ~ Age.baseline )
densityplot(mod_imp_qc, ~ Age )
densityplot(mod_imp_qc, ~ Age.fu )
densityplot(mod_imp_qc, ~ Educ )
densityplot(mod_imp_qc, ~ Ethn )
densityplot(mod_imp_qc, ~ WalPhyAct )
densityplot(mod_imp_qc, ~ ModPhyAct )
densityplot(mod_imp_qc, ~ VigPhyAct )
densityplot(mod_imp_qc, ~ PhyAct )
densityplot(mod_imp_qc, ~ PhyAct_bin )
densityplot(mod_imp_qc, ~ Smoke )
densityplot(mod_imp_qc, ~ Smoke_bin )
densityplot(mod_imp_qc, ~ Alcohol )
densityplot(mod_imp_qc, ~ Alcohol_bin )
densityplot(mod_imp_qc, ~ Diet_healthy )
densityplot(mod_imp_qc, ~ Diet_trad )
densityplot(mod_imp_qc, ~ Diet_proc.conf )
densityplot(mod_imp_qc, ~ Diet_proc)
densityplot(mod_imp_qc, ~ Diet_conf )
densityplot(mod_imp_qc, ~ Diet_semiveget )
densityplot(mod_imp_qc, ~ Diet_veget )
densityplot(mod_imp_qc, ~ Diet_healthy_bin )
densityplot(mod_imp_qc, ~ Diet_trad_bin )
densityplot(mod_imp_qc, ~ Diet_proc.conf_bin )
densityplot(mod_imp_qc, ~ Diet_semiveget_bin )
densityplot(mod_imp_qc, ~ Diet_proc_bin)
densityplot(mod_imp_qc, ~ Diet_conf_bin )
densityplot(mod_imp_qc, ~ Diet_veget_bin )

dev.off()

pdf(file = paste0('./results/convergence_imp-vs-obs_',date,'.pdf'))
plot(mod_imp_qc) 
dev.off()

# check for clear trends and whether chains do not mix well, i.e., there is more variation between the chains than within each chain
# for details and examples, see here: https://www.nerler.com/teaching/fgme2019/micheck


final_mod_imp = complete(mod_imp, action=20)
final_mod_imp = complete(mod_imp, action=30) # 30 for original imputation

col_order <- c("Age", "Sex", "Ethn",'Educ',
               "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
               'MDD','DM','CD','CMD','CVD1', 'PCM',
               'Alcohol','Alcohol_bin',
               'Smoke','Smoke_bin',
               'PhyAct','PhyAct_bin')

# fu ALSPAC partners
col_order <- c("Age", "Sex", "Ethn",'Educ',
               "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
               'MDD.baseline','MDD.fu',
               'DM.baseline','DM.fu',
               'CD',
               'CMD','CMD.fu',
               'CVD1',
               'PCM','PCM.fu',
               'Alcohol','Alcohol_bin',
               'Smoke','Smoke_bin',
               'PhyAct','PhyAct_bin')

# fu ALSPAC mothers
col_order <- c("Age.baseline","Age", "Sex", "Ethn",'Educ',
               "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
               'MDD.baseline','MDD.fu','MDD.fu1',
               'DM',
               'CD',
               'CMD',
               'CVD1',
               'PCM','PCM.fu',
               'Alcohol','Alcohol_bin',
               'Smoke','Smoke_bin',
               'PhyAct','PhyAct_bin')

# fu UKBB
col_order <- c("Age","Age.fu", "Sex", "Ethn",'Educ',
               "CM_EA", "CM_PA", 'CM_SA', 'CM','CM_sev',
               'MDD','MDD.fu',
               'DM_baseline','DM.fu',
               'CD','CD.fu',
               'CMD','CMD.fu',
               'CVD1','CVD1.fu',
               'PCM','PCM.fu',
               'Alcohol','Alcohol_bin',
               'Smoke','Smoke_bin',
               'PhyAct','PhyAct_bin')

final_mod_imp2 <- final_mod_imp[, c(col_order, grep("Diet_h", names(final_mod_imp), value=T),
                                    grep("veg", names(final_mod_imp), value=T),
                                    grep("Diet_con", names(final_mod_imp), value=T),
                                    grep("Diet_pro", names(final_mod_imp), value=T))]
names(final_mod_imp2)

final_mod_imp2 %>% 
  gtsummary::tbl_summary( 
    #by = condition,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),missing_text = "(Missing)", missing = 'always') %>%
  modify_caption("**Imputed Sample Characteristics**") %>%
  # this changes the table to a different format that can be saved as Word
  as_flex_table() %>% 
  # use flextable package to save table as word
  flextable::save_as_docx(path = paste0("./results/descriptives_Sample_imput",date,".docx"))#, pr_section = sect_properties)
