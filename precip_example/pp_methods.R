################################################################################
######### ISO CRPS decomposition for precipitation case study  #################
################################################################################

#-------------------------------------------------------------------------------
# load packages
library(Rcpp)
library(osqp)
library(devtools)
library(scoringRules)
library(ensembleMOS)
library(crch)
library(ensembleBMA)
library(isodistrreg)
#install_github("evwalz/isodisregSD")
library(isodisregSD)
#-------------------------------------------------------------------------------
# Set City: bru, lhr, zhr, fra

city <- 'bru'

# load data and functions to compute upper bound of the interval [0,b]
set_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
load(paste0(set_dir, "/precip_data/precipData_caseStudy.rda"))
data <- precipData_caseStudy
source(paste0(set_dir, "/functions_ab.R"))

for (lag in 1:5){
  #-------------------------------------------------------------------------------
  # read data for City
  data_sub <- subset(data, (airport==city)&(horizon==lag))
  n <- min(which(data_sub$date>="2015-01-01"))
  y <- data_sub$obs
  obs_train = data_sub$obs[1:(n-1)]
  obs_test = data_sub$obs[n:dim(data_sub)[1]]

  #-------------------------------------------------------------------------------
  # Postprocessing 1: HCLR

  grd <- seq(0, 150, 0.1)  # grid only used to make one CRPS computation!
  y <- data_sub$obs
  X <- data.frame( # transform everything to square root
    y = sqrt(data_sub$obs),
    hres = sqrt(data_sub$hres),
    ctr = sqrt(data_sub$ctr),
    mptb = apply(sqrt(data_sub[, paste0("p", 1:50)]), 1, mean),
    sd = apply(sqrt(data_sub[, paste0("p", 1:50)]), 1, sd))

  trDat = X[1:(n-1),,drop=FALSE]
  I = n:dim(X)[1]

  fit <- crch(
    formula = y ~ hres + ctr + mptb | sd,
    data = trDat,
    left = 0,
    type = "crps",
    dist = "logistic"
  )
  pred_HCLR <- predict(
    object = fit,
    newdata = X[I, , drop = FALSE],
    at = sqrt(grd),
    type = "probability"
  )
  paras_HCLR <- predict(
    object = fit,
    newdata = X[I, , drop = FALSE],
    #at = sqrt(grd),
    type = "parameter"
  )

  # Computation of the CRPS
  delta = grd[2]-grd[1]
  ind = matrix(NA,nrow=length(obs_test),ncol=length(grd))
  for (i in 1:length(obs_test)){
    ind[i,] <- obs_test[i]<=grd
  }

  crps_HCLR = mean(rowSums((pred_HCLR-ind)^2*delta))

  #computation of the upper bound of the interval [0,b] and the corresponding grid
  b= boundshclr_minmax(y=obs_test,loc=paras_HCLR[,1], scale=paras_HCLR[,2],
                       epsilon=crps_HCLR/1000,delta=diff(range(obs_test))/100)

  grid_hclr = seq(0,b,length.out=5000)


  # pred_HCLR on the grid [0,b]
  HCLR_cal_fit <- idrsd(y=obs_test, X = pred_HCLR, grid = grd, type = 'ecdf')


  #-------------------------------------------------------------------------------
  # Postprocessing 2: EMOS

  precipData_train = data_sub[1:(n-1),,drop=FALSE]
  precipData_test = data_sub[n:dim(data_sub)[1],,drop=FALSE]

  Data_train <- ensembleData(
    forecasts = precipData_train[,5:56],
    observations = precipData_train$obs,
    forecastHour = 24,
    #forecastHour = hori * 24L + 6L,
    initializationTime = 00,
    dates = format(precipData_train$date, "%Y%m%d"),
    exchangeable = setNames(
      c(1, 2, rep(3, 50)),
      c("hres", "ctr", paste0("p", 1:50))
    )
  )

  Data_test <- ensembleData(
    forecasts = precipData_test[,5:56],
    observations = precipData_test$obs,
    forecastHour = 24,
    #forecastHour = hori * 24L + 6L,
    initializationTime = 00,
    dates = format(precipData_test$date, "%Y%m%d"),
    exchangeable = setNames(
      c(1, 2, rep(3, 50)),
      c("hres", "ctr", paste0("p", 1:50))
    )
  )

  fit <- fitMOSgev0(Data_train)

  #  parameters: list of MEAN, SCALE and LOC
  Data_test_array <- as.matrix(Data_test[, 1:52])
  colnames(Data_test_array) <- NULL
  list_paras_emos <- paras_EMOS(fit,Data_test_array)

  crps_EMOS_ens = ensembleMOS::crps(fit=fit,ensembleData = Data_test)
  crps_EMOS <- mean(crps_EMOS_ens[, 2])

  b= bounds_pgev_minmax(obs_test,loc=list_paras_emos$LOC,scale=list_paras_emos$SCALE,
                        shape=list_paras_emos$SHAPE,epsilon=crps_EMOS/1000,delta=diff(range(obs_test))/100)

  grid_emos = seq(0,b,length.out=5000)

  pred_EMOS = ensembleMOS::cdf(fit=fit,ensembleData = Data_test,values=grid_emos)

  EMOS_cal_fit <- idrsd(y=obs_test, X = pred_EMOS, grid = grid_emos, type = 'ecdf')


  #-------------------------------------------------------------------------------
  # Postprocessing 3: BMA

  Data_train <- ensembleData(
    forecasts = cbind(
      hres = precipData_train$hres,
      ctr = precipData_train$ctr,
      mptb = apply(precipData_train[, paste0("p", 1:50)], 1, mean)
    ),
    observations = precipData_train$obs,
    forecastHour = 24,
    initializationTime = 00,
    dates = format(precipData_train$date, "%Y%m%d")
  )

  Data_test <- ensembleData(
    forecasts = cbind(
      hres = precipData_test$hres,
      ctr = precipData_test$ctr,
      mptb = apply(precipData_test[, paste0("p", 1:50)], 1, mean)
    ),
    observations = precipData_test$obs,
    forecastHour = 24,
    initializationTime = 00,
    dates = format(precipData_test$date, "%Y%m%d")
  )

  fit = fitBMAgamma0(Data_train)

  # parameters: list of MEAN, SCALE and LOC
  Data_test_array <- as.matrix(Data_test[, 1:3])
  colnames(Data_test_array) <- NULL
  list_paras_bma <- paras_BMA(fit,Data_test_array)

  crps_BMA_ens <- ensembleBMA::crps(fit=fit,ensembleData = Data_test)
  crps_BMA <- mean(crps_BMA_ens[, 2])

  b = bounds_bma_minmax(y=obs_test,weights=list_paras_bma$Weights,p0=list_paras_bma$P0,
                        mean=list_paras_bma$MEAN,variance=list_paras_bma$VAR,
                        epsilon=crps_BMA/1000,delta=diff(range(obs_test))/100)

  grid_bma= seq(0,b,length.out=5000)
  pred_BMA = ensembleBMA::cdf(fit=fit,ensembleData = Data_test,values=grid_bma)
  BMA_cal_fit <- idrsd(y=obs_test, X = pred_BMA, grid = grid_bma, type = 'ecdf')


  #-------------------------------------------------------------------------------
  # Postprocessing 4: IDR_cw

  X = data_sub[,5:56,drop=FALSE]
  X[,3]=apply(X[,3:52],1,mean)
  X = X[,1:3]
  colnames(X)= c("hres","ctr","ptm")
  data_train= X[1:(n-1),]
  data_test=X[n:dim(X)[1],]

  fit= isodistrreg::idr(obs_train,data_train)

  pred_IDR_cw <- isodistrreg:::predict.idrfit(fit, data = data_test)
  crps_IDR_cw = mean(isodistrreg:::crps(pred_IDR_cw, y = obs_test))

  IDR_cw_cal_fit <-idrsd(obs_test, X = pred_IDR_cw, type = 'idr')


  #-------------------------------------------------------------------------------
  # Postprocessing 5: IDR_st

  X = data_sub[,5:56,drop=FALSE]
  varNames = c("hres","ctr",paste0("p",1:50))
  groups = setNames(rep(1, 52), varNames)
  orders = c("sd" = 1)
  data_train= X[1:(n-1),]
  data_test=X[n:dim(X)[1],]

  fit = isodistrreg::idr(obs_train,data_train,orders=orders,groups=groups)

  pred_IDR_st <- isodistrreg:::predict.idrfit(fit, data = data_test,orders=orders,groups=groups)
  crps_IDR_st = mean(isodistrreg:::crps(pred_IDR_st, y = obs_test))

  IDR_st_cal_fit <- idrsd(obs_test, X = pred_IDR_st, type = 'idr')

  #-------------------------------------------------------------------------------
  # CPRS of the post-processed forecasts
  crps_PP = c(crps_HCLR, crps_EMOS,crps_BMA,crps_IDR_cw,crps_IDR_st)

  #-------------------------------------------------------------------------------
  # CRPS of the calibrated post-processed forecasts

  crps_PP_cal = c(
    mean(isodisregSD:::crps(predict(HCLR_cal_fit), y = obs_test)),
    mean(isodisregSD:::crps(predict(EMOS_cal_fit), y = obs_test)),
    mean(isodisregSD:::crps(predict(BMA_cal_fit), y = obs_test)),
    mean(isodisregSD:::crps(predict(IDR_cw_cal_fit), y = obs_test)),
    mean(isodisregSD:::crps(predict(IDR_st_cal_fit), y = obs_test)))

  # CPRS of the reference forecast

  crps_ref = rep(mean(crps_sample(obs_test,matrix(rep(obs_test,length(obs_test)),byrow=TRUE,nrow =length(obs_test)))),5)

  #-------------------------------------------------------------------------------
  # Score Decomposition

  mcb = crps_PP - crps_PP_cal
  dsc = crps_ref - crps_PP_cal
  unc = crps_ref
  decomposition = data.frame(rbind(mcb,dsc,unc))
  crps_values = data.frame(rbind(crps_PP,crps_PP_cal,crps_ref))
  colnames(decomposition) =
    colnames(crps_values) = c("HCLR","EMOS","BMA","IDR_cw","IDR_st")

  results = data.frame(matrix(NA,ncol=5,nrow=20))
  colnames(results)=c("airport","horizon","method","category","value")
  results[,1]= city
  results[,2]= lag
  results$value[1:5]=crps_PP
  results$category[1:5]="crps"
  results$value[6:10]=mcb
  results$category[6:10]="mcb"
  results$value[11:15]=dsc
  results$category[11:15]="dsc"
  results$value[16:20]=unc
  results$category[16:20]="unc"
  results$method = rep(c("HCLR","EMOS","BMA","IDR_cw","IDR_st"),4)


  export=list("decomposition"=decomposition,"crps_values"=crps_values,
              "data"=results)
  save(export,file=paste0(set_dir, "/deco_data/pp_",city,"_lag=", lag,".rda"))
}
