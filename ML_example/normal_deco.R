# Methods with predictive normal distributions:

# 1. Single gaussian: Same variance accross samples, thus use "sd" on predictive normal distributions
# 2. Laplace: Different variance, thus use "sd" on predictive normal distributions on Interval [a, b]


#############################################################################################
#############################################################################################

library(isodisregSD)

# Load function script
source('functions_deco.R')


# Loop over methods and data sets
methods <- c( 'laplace', 'single_gaussian')
data_dirs <- c('bostonHousing', 'concrete', 'energy', 'yacht')

for (method in methods){
  y_file = 'y_test_split_'
  if (method == 'single_gaussian'){
    mean_file <- 'fct_test_10_xepochs_1_hidden_layers_'
    y_file <- 'y_test_10_xepochs_1_hidden_layers_'
  } else {
    mean_file <- 'fct_test_norm_10_xepochs_1_hidden_layers_'
    std_file <- 'fct_test_sg_10_xepochs_1_hidden_layers_'
    y_file <- 'y_test_10_xepochs_1_hidden_layers_'
  }
  

  for (data_directory in data_dirs){ # , 'power-plant'
    deco_path <- paste('./deco_data/', method, '/UCI_datasets/',  data_directory , '/deco/', sep = '')
    
    dim_split <- 20
    if (data_directory == 'protein-tertiary-structure'){
      dim_split <- 5
    }
    
    eps <- "sd"
    CRPS <- rep(0, dim_split)
    DSC <- rep(0, dim_split)
    MSC <- rep(0, dim_split)
    UNC <- rep(0, dim_split)
    
    if (method == 'single_gaussian'){
      if (data_directory == 'protein-tertiary-structure') {
        sigma_all <- read.table(paste(deco_path, 'sigma_test_10_xepochs_1_hidden_layers_4_split.txt', sep = ''))
      } else {
        sigma_all <- read.table(paste(deco_path, 'sigma_test_10_xepochs_1_hidden_layers_19_split.txt', sep = ''))
      }
      col_string <- colnames(sigma_all)
      sigma_all <- sigma_all[, col_string]
    }
    
    for (i in 1:dim_split){
      n_split <- i-1
      if (method == 'single_gaussian'){
        y_test <-  read.table(paste(deco_path,y_file, n_split,'_split.txt', sep = ''))
        mean <- read.table(paste(deco_path,mean_file, n_split,'_split.txt', sep = ''))
        col_string <- colnames(y_test)
        y_test <- y_test[, col_string]
        mean <- mean[, col_string]
        sigma <- rep(sigma_all[i], length(y_test))
      } else {
        y_test <-  read.table(paste(deco_path,y_file, n_split,'_split.txt', sep = ''))
        mean <- read.table(paste(deco_path,mean_file, n_split,'_split.txt', sep = ''))
        sigma <- read.table(paste(deco_path,std_file, n_split,'_split.txt', sep = ''))
        col_string <- colnames(y_test)
        y_test <- y_test[, col_string]
        mean <- mean[, col_string]
        sigma <- sigma[, col_string]
      }

      paras <- matrix(0, nrow = length(y_test), ncol = 2)
      paras[, 1] <- mean
      paras[, 2] <- sigma

      crps_gauss <- gaussian_crps(y_test, mean, sigma)
      
      if (method == 'single_gaussian'){
        cali_idr <-  idrsd(y_test, X = paras, type = 'normal')
        cali_preds <- predict(cali_idr)
        cali_crps <- mean(crps(cali_preds, y_test))
        
        uncertainty <- crps_unc(y_test)
        
        MSC[i] <- crps_gauss - cali_crps
        DSC[i] <- uncertainty - cali_crps
        UNC[i] <- uncertainty
        CRPS[i] <- crps_gauss

      } else {
        delta <- (max(y_test) - min(y_test)) / 100
        epsilon <- 0.001 * crps_gauss
        
        interval <- bounds_norm_minmax(y_test, mean, sigma, epsilon, delta)
        a = interval[1]
        b <- interval[2]
        
        cali_idr <-  idrsd(y_test, X = paras, eps ='sd', type = 'normal_ab', inta = a, intb = b)
        cali_preds <- predict(cali_idr)
        cali_crps <- mean(crps(cali_preds, y_test))
        
        uncertainty <- crps_unc(y_test)
        
        crps_norm_ab <- func_crps_normAB(y_test, mean, sigma, a, b)
        
        MSC[i] <- crps_norm_ab - cali_crps
        DSC[i] <- uncertainty - cali_crps
        UNC[i] <- uncertainty
        CRPS[i] <- crps_norm_ab
      }
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
    colnames(mean_scores) <- NULL
    write.table(mean_scores, paste(deco_results_path, "log_scores.txt", sep=""), col.names = FALSE)
  }
}



