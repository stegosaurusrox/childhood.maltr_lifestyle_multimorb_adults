# CREDIT:
# some of this code was written by
# Serena Defina: https://github.com/SereDef
# and adjusted, merged, extended by Esther Walton

library(openxlsx)
library(mice)

date <- format(Sys.Date(), "%d%m%y")


# SAMPLE (AFTER mod_impATION ) =================================

describe <- function(mod_imp=mod_imp) {
  # determine categorical and continuous vars 
  lvl_length <- lapply(mod_imp$data, function(var) length(levels(as.factor(var)))) 
  # Cutoff 15 levels: consider it categorical
  cat_vars <- names(which(lvl_length < 5))
  
  # Stack mod_imped datasets in long format, excluding the original data
  impdat <- mice::complete(mod_imp, action="long", include = F)
  # Set to factors or numeric when appropiate
  impdat[,cat_vars] <- lapply(impdat[,cat_vars] , as.factor)
  impdat = impdat[,-grep('ID',names(impdat))] # ID for UKK, cidB2957 for ALSPAC
  impdat[, !names(impdat)%in%cat_vars] <- lapply(impdat[,!names(impdat)%in%cat_vars] , as.numeric)
  
  pool_descriptives <- function(implist, column_names, categorical=T) {
    summ <- with(implist, by(implist, .imp, function(x) summary(x[, -c(1,2)],digits=4))) 
    if (categorical==F) {
      # Pool summary 
      num_pool <- lapply(summ, function(m) matrix(as.numeric(sapply(strsplit(m, ":"), "[[", 2)), nrow = dim(m)[1], ncol=dim(m)[2]))
      pool_mean <- Reduce("+",num_pool)/length(num_pool)
      # Pool SDs
      sds <- with(implist, by(implist, .imp, function(x) round(apply(x[, -c(1, 2)], 2, sd, na.rm = T), 4)))
      pool_sds <- Reduce("+",sds)/length(sds)
      # Bind SDs to other metrics
      summ_df <- data.frame(rbind(pool_mean,pool_sds))
      # Define column and row names
      colnames(summ_df) <- colnames(implist[-c(1,2)])
      rownames(summ_df) <- c('Min','1stQ','Median','Mean','3rdQ','Max','SD')
    } else { 
      pool_mean <- Reduce("+",summ)/length(summ) 
      summ_df <- data.frame(pool_mean)
      colnames(summ_df) <- 'counts' # colnames(implist[-c(1,2)])
      rownames(summ_df) <- names(pool_mean)
    }
    return(summ_df)
  }
  # Continuous data
  cnt <- impdat[, -c(which(colnames(impdat) %in% cat_vars))]
  cnt_summ <- cbind(c('Min','1stQ','Median','Mean','3rdQ','Max','SD'), 
                    pool_descriptives(cnt, categorical = F))
  
  # Categorical / binary
  cat_summ <- NA
  for (v in cat_vars) {
    sel <- impdat[, c('.imp','.id', v)]
    v_summ <- pool_descriptives(sel)
    v_summ <- cbind(row.names(v_summ), v_summ)
    v_summ$percent <- (v_summ$count / nrow(mod_imp$data))*100
    cat_summ <- rbind(cat_summ, v, v_summ, NA)
  }
  
  # Correlation matrix in the mod_imped set
  cors_imp <- miceadds::micombine.cor(mi.res = impdat, method='spearman',
                                      variables = colnames(impdat)[!colnames(impdat) %in% 
                                               c('.imp','.id',cat_vars)]) 
  
  # Export the outputs of summary statistics into an xlsx file with one model per sheet
  stats <- list('s_imp_cnt' = cnt_summ, 's_imp_cat' = cat_summ, 'cor_imp' = cors_imp)
  
  return(stats)
}

s = describe(mod_imp)


openxlsx::write.xlsx(s,file = paste0('./results/ALSPAC_Descriptives_pooled_',date,'.xlsx'),overwrite=T)
