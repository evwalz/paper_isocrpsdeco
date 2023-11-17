################################################################################
######### ISO CRPS decomposition for precipitation case study  #################
################################################################################

#-------------------------------------------------------------------------------
# load packages
library(Rcpp)
library(osqp)
library(devtools)
#install_github("evwalz/isodisregSD")
library(isodisregSD)

city <- 'bru' # bru, lhr, zhr, fra

set_dir = dirname(rstudioapi::getSourceEditorContext()$path)
load(paste0(set_dir, "/data/precipData_caseStudy.rda"))
data <- precipData_caseStudy

for (lag in 1:5){
  #-------------------------------------------------------------------------------
  # read data for City
  data_sub <- subset(data, (airport==city)&(horizon==lag))
  n <- min(which(data_sub$date>="2015-01-01"))
  y <- data_sub$obs
  obs_train = data_sub$obs[1:(n-1)]
  obs_test = data_sub$obs[n:dim(data_sub)[1]]
  X = data_sub[,5:56,drop=FALSE]
  data_train= X[1:(n-1),]
  data_test=X[n:dim(X)[1],]
  #-------------------------------------------------------------------------------
  # Decomposition for ensemble forecasts

  deco = isodeco_crps(y=obs_test, X = data_test,type="ensemble")
  mcb = deco$MCB
  dsc = deco$DSC
  unc = crps_ref=deco$UNC
  crps_ens = deco$CRPS
  crps_ens_cal=crps_ens-mcb
  #-------------------------------------------------------------------------------
  # Generate Dataframe to export

  decomposition = data.frame(cbind(mcb,dsc,unc))
  crps_values = data.frame(cbind(crps_ens,crps_ens_cal,crps_ref))

  results = data.frame(matrix(NA,ncol=5,nrow=4))
  colnames(results)=c("airport","horizon","method","category","value")
  results[,1]= city
  results[,2]= lag
  results$value[1]=deco$CRPS
  results$category[1]="crps"
  results$value[2]=mcb
  results$category[2]="mcb"
  results$value[3]=dsc
  results$category[3]="dsc"
  results$value[4]=unc
  results$category[4]="unc"
  results$method = "ens"

  export=list("decomposition"=decomposition,"crps_values"=crps_values,
              "data"=results)

  save(export,file=paste0(set_dir, "/deco_data/ens_",city,"_lag=", lag,".rda"))
}
