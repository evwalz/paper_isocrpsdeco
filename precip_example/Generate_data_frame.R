set_dir = dirname(rstudioapi::getSourceEditorContext()$path)

Data= data.frame(matrix(NA,ncol=5,nrow=1))
colnames(Data)=c("airport","horizon","method","category","value")
for (i in 1:5){
  load(paste0(set_dir, "/deco_data/pp/pp_lhr_lag=",i,".rda"))
  Data=rbind(Data,export$data)
  load(paste0(set_dir, "/deco_data/ens/ens_lhr_lag=",i,".rda"))
  Data=rbind(Data,export$data)

}
for (i in 1:5){
    load(paste0(set_dir, "/deco_data/pp/pp_zrh_lag=",i,".rda"))
    Data=rbind(Data,export$data)
    load(paste0(set_dir, "/deco_data/ens/ens_zrh_lag=",i,".rda"))
    Data=rbind(Data,export$data)
}
for (i in 1:5){
    load(paste0(set_dir, "/deco_data/pp/pp_bru_lag=",i,".rda"))
    Data=rbind(Data,export$data)
    load(paste0(set_dir, "/deco_data/ens/ens_bru_lag=",i,".rda"))
    Data=rbind(Data,export$data)
}

for (i in 1:5){
    load(paste0(set_dir, "/deco_data/pp/pp_fra_lag=",i,".rda"))
    Data=rbind(Data,export$data)
    load(paste0(set_dir, "/deco_data/ens/ens_fra_lag=",i,".rda"))
    Data=rbind(Data,export$data)
}
Data= Data[-1,]
save(Data,file=paste0(set_dir, "/deco_data/results.rda"))
