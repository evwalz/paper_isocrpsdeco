################################################################################
############## Comparing PP Methods with CRPS Decomposition ####################
################################################################################

#-------------------------------------------------------------------------------
# load packages
library(Rcpp)
library(osqp)

set_dir = dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(set_dir, '/functions_deco.R', sep = ''))

#-------------------------------------------------------------------------------
# Set city and lag
city <- 'bru' # 'bru', 'lhr', 'zrh', 'fra'
lag <- 1 # 1, 2, 3, 4, 5

#-------------------------------------------------------------------------------
# read data
load("/Users/eva/Documents/Work/promotion/decomposition/docs/crps_decomposition/Old_files/Case Study/precipData_caseStudy.rda")
load("precipData_caseStudy.rda")
data <- precipData_caseStudy
data_sub <- subset(data, (airport==city)&(horizon==lag))
n <- min(which(data_sub$date>="2015-01-01"))
y <- data_sub$obs

obs_train = data_sub$obs[1:(n-1)]
obs_test = data_sub$obs[n:dim(data)[1]]


#-------------------------------------------------------------------------------
# Postprocessing 1: HCLR


#-------------------------------------------------------------------------------
# Postprocessing 2: EMOS


#-------------------------------------------------------------------------------
# Postprocessing 3: BMA


#-------------------------------------------------------------------------------
# Postprocessing 4: IDR_cw

#-------------------------------------------------------------------------------
# Postprocessing 5: IDR_st


#-------------------------------------------------------------------------------
# Calculation of the CPRS of the original post-processed forecasts


#-------------------------------------------------------------------------------
# IDRsd on HCLR, EMOS, BMA, IDR_cw and IDR_st


#-------------------------------------------------------------------------------
# Calculation of the CPRS of the recalibrated post-processed forecasts


#-------------------------------------------------------------------------------
# Score Decomposition

mcb = crps_PP - crps_PP_cal
dsc = crps_ref - crps_PP_cal
unc = crps_ref
decomposition = data.frame(rbind(mcb,dsc,unc))
crps_values = data.frame(rbind(crps_PP,crps_PP_cal,crps_ref))
colnames(decomposition) =
  colnames(crps_values) = c("HCLR","EMOS","BMA","IDR_cw","IDR_afsd")

data = data.frame(matrix(NA,ncol=5,nrow=20))
colnames(data)=c("airport","horizon","method","category","value")
data[,1]= city
data[,2]= lag
data$value[1:5]=crps_PP
data$category[1:5]="crps"
data$value[6:10]=mcb
data$category[6:10]="mcb"
data$value[11:15]=dsc
data$category[11:15]="dsc"
data$value[16:20]=unc
data$category[16:20]="unc"
data$method = rep(c("HCLR","EMOS","BMA","IDR_cw","IDR_st"),4)


export=list("decomposition"=decomposition,"crps_values"=crps_values,
            "data"=data)

save(export,file=paste0("pp_bru=",city,"_lag=", lag ,".rda"))

