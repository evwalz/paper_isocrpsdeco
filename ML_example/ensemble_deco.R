# all ensemble/ecdf methods:

# 1. EasyUQ
# 2. CP

library(extraDistr)
library(scoringRules)

source('functions_deco.R')
######################################################################

method = 'easyuq' 
epochs_multiplier <- 10
num_hidden_layers <- 1

data_dirs <- c('bostonHousing', 'energy', 'yacht', 'concrete') # ,'kin8nm'

for (data_directory in data_dirs){ # , 'power-plant'
  deco_path <- paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/deco/', sep = '')
  
  dim_split <- 20
  if (data_directory == 'protein-tertiary-structure'){
    dim_split <- 5
  }
  
  CRPS <- rep(0, dim_split)
  DSC <- rep(0, dim_split)
  MSC <- rep(0, dim_split)
  UNC <- rep(0, dim_split)
  
  
  for (i in 1:dim_split){
    n_split <- i-1
    
    y_test <- read.table(paste(deco_path,'y_test_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    y_train <- read.table(paste(deco_path,'y_train_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    fct_test <- read.table(paste(deco_path,'fct_test_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    fct_train <- read.table(paste(deco_path,'fct_train_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    
    idr_train <- isodistrreg::idr(y = y_train, X = data.frame('fct' = fct_train))
    idr_preds <-predict(idr_train, data.frame('fct' = fct_test))
    crps_val <- mean(crps(idr_preds, y_test))
    
    idr_cali <- idrsd(y = y_test, X = idr_preds,type = 'idr')
    cali_crps <- mean(crps2(predict(idr_cali), y = y_test))
    uncertainty <- crps_unc(y_test)
    
    MSC[i] <- crps_val - cali_crps
    DSC[i] <- uncertainty - cali_crps
    UNC[i] <- uncertainty
    CRPS[i] <- crps_val
    
  }
  
  nMSC <- 1 - (1 / UNC)*MSC
  nDSC <- (1 / UNC)*DSC
  deco_results_path <- paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/deco_results/', sep = '')
  
  if (file.exists(deco_results_path) == FALSE) {
    dir.create(deco_results_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  write.table(MSC, paste(deco_results_path, "msc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(DSC, paste(deco_results_path, "dsc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(UNC, paste(deco_results_path, "unc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(CRPS, paste(deco_results_path, "crps.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(nMSC, paste(deco_results_path, "norm_msc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(nDSC, paste(deco_results_path, "norm_dsc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  
  mean_scores <- data.frame(c(mean(MSC), mean(DSC), mean(UNC), mean(CRPS), mean(nMSC), mean(nDSC)))
  
  row.names(mean_scores) <- c('MSC ', 'DSC ', 'UNC ', 'CRPS', 'nMSC', 'nDSC')
  write.table(mean_scores, paste(deco_results_path, "log_scores.txt", sep=""), col.names = FALSE)
}




##################################### CP ########################################

method = 'cp'


ens_file <- 'ens_test_10_xepochs_1_hidden_layers_'
y_file <- 'y_test_10_xepochs_1_hidden_layers_'

data_dirs <- c('bostonHousing', 'concrete', 'energy', 'yacht')
for (data_directory in data_dirs){ 
  deco_path <- paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/deco/', sep = '')
  
  dim_split <- 20
  # overall splits
  if (data_directory == 'protein-tertiary-structure'){
    dim_split <- 5
  }
  
  CRPS <- rep(0, dim_split)
  DSC <- rep(0, dim_split)
  MSC <- rep(0, dim_split)
  UNC <- rep(0, dim_split)

  for (i in 1:dim_split){
    n_split <- i-1
    
    y_test <-  read.table(paste(deco_path,y_file, n_split,'_split.txt', sep = ''))
    ens <- read.table(paste(deco_path,ens_file, n_split,'_split.txt', sep = ''))
    col_string <- colnames(y_test)
    y_test <- y_test[, col_string]
    
    cali_idr <-  idr2(y_test, ens, type = 'ensemble')
    cali_preds <- predict(cali_idr)
    
    cali_crps <- mean(crps(cali_preds, y_test))
    
    uncertainty <- crps_unc(y_test)
    ens_mat <- as.matrix(ens)
    colnames(ens_mat) <- NULL
    rownames(ens_mat) <- NULL
    crps_val <- mean(crps_sample(y_test, ens_mat))
    
    MSC[i] <- crps_val - cali_crps
    DSC[i] <- uncertainty - cali_crps
    UNC[i] <- uncertainty
    CRPS[i] <- crps_val
  }
  
  # save values to folder
  nMSC <- 1 - (1 / UNC)*MSC
  nDSC <- (1 / UNC)*DSC
  deco_results_path <- paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/deco_results/', sep = '')
  
  if (file.exists(deco_results_path) == FALSE) {
    dir.create(deco_results_path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  write.table(MSC, paste(deco_results_path, "msc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(DSC, paste(deco_results_path, "dsc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(UNC, paste(deco_results_path, "unc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(CRPS, paste(deco_results_path, "crps.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(nMSC, paste(deco_results_path, "norm_msc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  write.table(nDSC, paste(deco_results_path, "norm_dsc.txt", sep=""), row.names=FALSE, col.names = FALSE)
  
  mean_scores <- data.frame(c(mean(MSC), mean(DSC), mean(UNC), mean(CRPS), mean(nMSC), mean(nDSC)))
  
  row.names(mean_scores) <- c('MSC ', 'DSC ', 'UNC ', 'CRPS', 'nMSC', 'nDSC')
  write.table(mean_scores, paste(deco_results_path, "log_scores.txt", sep=""), col.names = FALSE)
}


