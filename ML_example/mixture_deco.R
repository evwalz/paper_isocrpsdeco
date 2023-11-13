# all mixture distribution methods:

# 1. smooth EasyUQ
# 2. smooth CP
# 3. MC Dropout

library(extraDistr)

source('functions_deco.R')

epochs_multiplier <- 10
num_hidden_layers <- 1


method = 'easyuq_smooth' 
data_dirs <- c('bostonHousing', 'concrete', 'energy', 'yacht')


for (data_directory in data_dirs){
  deco_path <- paste('./deco_data/', substring(method, 1, 6), '/UCI_datasets/',  data_directory , '/deco/', sep = '')
  smooth_path <- paste('./deco_data/', substring(method, 1, 6), '/UCI_datasets/',  data_directory , '/results/', sep = '')
  
  dim_split <- 20
  if (data_directory == 'protein-tertiary-structure'){
    dim_split <- 5
  }
  
  CRPS <- rep(0, dim_split)
  DSC <- rep(0, dim_split)
  MSC <- rep(0, dim_split)
  UNC <- rep(0, dim_split)
  
  # CRPS of smoothed object
  df_all <- read.table(paste(smooth_path,'test_df_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers.txt', sep = ''))$V1
  h_all <- read.table(paste(smooth_path,'test_bw_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers.txt', sep = ''))$V1
  
  #i <- 1
  for (i in 1:dim_split){
    n_split <- i-1
    
    y_test <- read.table(paste(deco_path,'y_test_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    y_train <- read.table(paste(deco_path,'y_train_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    fct_test <- read.table(paste(deco_path,'fct_test_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    fct_train <- read.table(paste(deco_path,'fct_train_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    
    idr_train <- isodistrreg::idr(y = y_train, X = data.frame('fct' = fct_train))
    idr_preds <-predict(idr_train, data.frame('fct' = fct_test))
    df <- df_all[i]
    h <- h_all[i]
    
    delta <- (max(y_test) - min(y_test)) / 100
    
    crps_original <- crps(idr_preds, y_test)
    
    if (df == "None"){
      crps_val = rep(0, length(y_test))
      epsilon <- 0.001 * crps_original
      interval <- bounds_norm_mix(y_test,idr_preds, h, epsilon, delta)
      grid_vals <- seq(interval[1], interval[2], length.out = 5000)

      pmix_ecdf_test <- matrix(, nrow =  length(y_test), ncol = length(grid_vals))
      for (j in 1:length(y_test)){
        mean <- idr_preds[[j]]$points
        weights <- diff(c(0, idr_preds[[j]]$cdf))
        for (k in 1:length(grid_vals)){
          pmix_ecdf_test[j, k] <- extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
        }
      }
    } else {
      df = as.numeric(df)
      crps_val = rep(0, length(y_test))
      epsilon <- 0.001 * crps_original

      interval <- bounds_t_mix(y_test,idr_preds, h,df, epsilon, delta)
      grid_vals <- seq(interval[1], interval[2], length.out = 5000)
      pmix_ecdf_test <- matrix(, nrow =  length(y_test), ncol = length(grid_vals))
      for (j in 1:length(y_test)){
        mean <- idr_preds[[j]]$points
        weights <- diff(c(0, idr_preds[[j]]$cdf))
        for (k in 1:length(grid_vals)){
          pmix_ecdf_test[j, k] <- sum(weights * pt((grid_vals[k] - mean) / h, df = df))
        }
      }
    }
    
    crps_val <- crps_ecdf(y_test, pmix_ecdf_test, grid_vals)
    
    smooth_cali <- idrsd(y_test, X = pmix_ecdf_test, grid = grid_vals, type = 'ecdf')
    #idr_cali <- idr2(y = y_test, X = ,type = 'dis', eps = 'sd')
    cali_crps <- mean(crps.idr(predict(smooth_cali), y = y_test))
    
    uncertainty <- crps_unc(y_test)
    
    MSC[i] <- mean(crps_val) - cali_crps
    DSC[i] <- uncertainty - cali_crps
    UNC[i] <- uncertainty
    CRPS[i] <- mean(crps_val)
    
  }
  
  # save values to folder
  nMSC <- 1 - (1 / UNC)*MSC
  nDSC <- (1 / UNC)*DSC
  deco_results_path <- paste('./deco_data/', substring(method, 1, 6), '/UCI_datasets/',  data_directory, '/deco_results_smooth/', sep = '')
  
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



###################################### CP #################################################

method = 'cp_smooth'


ens_file <- 'ens_test_10_xepochs_1_hidden_layers_'
y_file <- 'y_test_10_xepochs_1_hidden_layers_'

data_dirs <- c('bostonHousing','energy', 'concrete', 'yacht')
for (data_directory in data_dirs){ 
  deco_path <- paste('./deco_data/', substring(method, 1, 2), '/UCI_datasets/',  data_directory , '/deco/', sep = '')
  deco_path2 <- paste('./deco_data/', substring(method, 1, 2), '/UCI_datasets/',  data_directory , '/results/', sep = '')
  
  dim_split <- 20
  # overall splits
  if (data_directory == 'protein-tertiary-structure'){
    dim_split <- 5
  }
  
  CRPS <- rep(0, dim_split)
  DSC <- rep(0, dim_split)
  MSC <- rep(0, dim_split)
  UNC <- rep(0, dim_split)
  
  df_all <- read.table(paste(deco_path2, 'test_df_10_xepochs_1_hidden_layers.txt', sep =''))$V1
  h_all <- read.table(paste(deco_path2, 'test_bw_10_xepochs_1_hidden_layers.txt', sep=''))$V1
  
  crps_all <- read.table(paste(deco_path2,'test_crps_smooth_10_xepochs_1_hidden_layers.txt', sep = ''))$V1
  
  for (i in 1:dim_split){
    n_split <- i-1
    
    y_test <-  read.table(paste(deco_path,y_file, n_split,'_split.txt', sep = ''))
    ens <- read.table(paste(deco_path,ens_file, n_split,'_split.txt', sep = ''))
    ens_mat <- as.matrix(ens)
    rownames(ens_mat) <- NULL
    colnames(ens_mat) <- NULL
    col_string <- colnames(y_test)
    y_test <- y_test[, col_string]
    
    df <- df_all[i]
    h <- h_all[i]
    delta <- (max(y_test) - min(y_test)) / 100
    
    crps_original <- crps_all[i]
    
    if (df == "None"){
      crps_val = rep(0, length(y_test))
      epsilon <- 0.001 * crps_original
      
      interval <- bounds_norm_mix_cp(y_test,ens_mat, h, epsilon, delta)
      
      grid_vals <- seq(interval[1], interval[2], length.out = 5000)
      #grid_vals <- seq(0, 100, 0.1)
      pmix_ecdf_test <- matrix(, nrow =  length(y_test), ncol = length(grid_vals))
      for (j in 1:length(y_test)){
        mean <- sort(unique(ens[j,]))
        colnames(mean) <- NULL
        mean <- unlist(mean)
        ecdf_fun <- ecdf(ens_mat[j,])
        ecdf_vals <- ecdf_fun(mean)
        weights <-  diff(c(0, ecdf_vals))
        
        for (k in 1:length(grid_vals)){
          pmix_ecdf_test[j, k] <- extraDistr::pmixnorm(grid_vals[k], mean = mean, alpha = weights, sd = rep(h, length(mean)))
        }
      }
    } else {
      df = as.numeric(df)
      crps_val = rep(0, length(y_test))
      epsilon <- 0.001 * crps_original
      
      interval <- bounds_t_mix_cp(y_test,ens_mat, h, df, epsilon, delta)
      grid_vals <- seq(interval[1], interval[2], length.out = 5000)
      pmix_ecdf_test <- matrix(, nrow =  length(y_test), ncol = length(grid_vals))
      for (j in 1:length(y_test)){
        mean <- sort(unique(ens[j,]))
        colnames(mean) <- NULL
        mean <- unlist(mean)
        ecdf_fun <- ecdf(ens_mat[j,])
        ecdf_vals <- ecdf_fun(mean)
        weights <-  diff(c(0, ecdf_vals))
        for (k in 1:length(grid_vals)){
          pmix_ecdf_test[j, k] <- sum(weights * pt((grid_vals[k] - mean) / h, df = df))
        }
      }
    }

    crps_val <- crps_ecdf(y_test, pmix_ecdf_test, grid_vals)
    
    smooth_cali <- idrsd(y_test, X = pmix_ecdf_test, grid = grid_vals, type = 'ecdf')
    cali_crps <- mean(crps.idr(predict(smooth_cali), y = y_test))
    
    uncertainty <- crps_unc(y_test)
    
    MSC[i] <- mean(crps_val) - cali_crps
    DSC[i] <- uncertainty - cali_crps
    UNC[i] <- uncertainty
    CRPS[i] <- mean(crps_val)
    
  }
  
  
  # save values to folder
  nMSC <- 1 - (1 / UNC)*MSC
  nDSC <- (1 / UNC)*DSC
  deco_results_path <- paste('./deco_data/', substring(method, 1, 2), '/UCI_datasets/',  data_directory, '/deco_results_smooth/', sep = '')
    
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




###################################### MC #################################################

epochs_multiplier <- 10
num_hidden_layers <- 1

method = 'mc_dropout'

data_dirs <- c('bostonHousing', 'energy', 'concrete', 'yacht')
for (data_directory in data_dirs){ 
  deco_path <- paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/deco/', sep = '')
  tau_path <-paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/results/', sep = '')
  
  dim_split <- 20
  if (data_directory == 'protein-tertiary-structure'){
    dim_split <- 5
  }
  
  CRPS <- rep(0, dim_split)
  DSC <- rep(0, dim_split)
  MSC <- rep(0, dim_split)
  UNC <- rep(0, dim_split)
  
  tau_vals <- read.table(paste(tau_path,'test_tau_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers.txt', sep = ''))$V1
  
  for (i in 1:dim_split){
    n_split <- i-1
    tau = sqrt(tau_vals[i])
    y_test <- read.table(paste(deco_path,'y_test_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = ''))$V1
    y_that <- as.matrix(read.table(paste(deco_path,'Ythat_test_',epochs_multiplier ,'_xepochs_',num_hidden_layers ,'_hidden_layers_', n_split,'_split.txt', sep = '')))
    colnames(y_that) <- NULL
    y_that2 <- t(y_that)
    sigmas <- matrix( 1 /tau, ncol = dim(y_that2)[2], nrow = dim(y_that2)[1])
    
    crps_mc <- mean(scoringRules::crps_mixnorm(y_test, m = y_that2, s =sigmas))
    
    
    delta <- (max(y_test) - min(y_test)) / 100
    epsilon <- 0.001*crps_mc
    
    
    crps_val = rep(0, length(y_test))
    
    interval <- bounds_norm_mix_mc(y_test,y_that2, sigmas, epsilon, delta)
    grid_vals <- seq(interval[1], interval[2], length.out = 5000)
    
    pmix_ecdf_test <- matrix(, nrow =  length(y_test), ncol = length(grid_vals))
    for (j in 1:length(y_test)){
      for (k in 1:length(grid_vals)){
        pmix_ecdf_test[j, k] <- extraDistr::pmixnorm(grid_vals[k], mean = y_that2[j,], sd = sigmas[j,], alpha = rep(1, length(sigmas[j,]))/ length(sigmas[j,]))
      }
    }
    
    crps_val <- crps_ecdf(y_test, pmix_ecdf_test, grid_vals)
    
    smooth_cali <- idr2(y_test, X = pmix_ecdf_test, grid = grid_vals, eps = 'sd', type = 'ecdf')
    
    cali_crps <- mean(crps.idr(predict(smooth_cali), y = y_test))
    
    uncertainty <- crps_unc(y_test)
    
    print( mean(crps_val) - cali_crps)
    
    MSC[i] <- mean(crps_val) - cali_crps
    DSC[i] <- uncertainty - cali_crps
    UNC[i] <- uncertainty
    CRPS[i] <- mean(crps_val)
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
