


normalizer_ref= function(dat, ref_name){
  microbiome=dat %>% select(-c(Sample.ID,pop))
  nbig=dim(microbiome)[1];    d=dim(microbiome)[2]
  datasplit<-list();    sid_list=list()
  for (i in 1:length(unique(dat$pop))){ datasplit[[i]]=filter(dat, pop==unique(dat$pop)[i])
  sid_list[[i]]=datasplit[[i]]$Sample.ID;      datasplit[[i]]=datasplit[[i]] %>% select(-c(Sample.ID,pop))
  datasplit[[i]]=datasplit[[i]]+1};
  mean_ref_logscale=matrix(0,length(unique(dat$pop)),1);    normalize_logscale=list() 
  for (i in 1:length(unique(dat$pop))){
    mean_ref_logscale[i]=mean(log(as.matrix(datasplit[[i]][,ref_name])),na.rm=T)
    normalize_logscale[[i]]=log(as.matrix(datasplit[[i]][,ref_name]))-mean_ref_logscale[i]}
  normalize_ref=cbind(sid_list[[1]],normalize_logscale[[1]])
  for (i in 2:length(unique(dat$pop))){
    temp=cbind(sid_list[[i]],normalize_logscale[[i]]);      normalize_ref=rbind(normalize_ref,temp)}
  normalize_ref=data.frame(normalize_ref);   names(normalize_ref)=c("Sample.ID","normal")
  return(normalize_ref)}
  
  



  normalizer_gm=function(dat){
    microbiome=dat %>% select(-c(Sample.ID,pop))
    nbig=dim(microbiome)[1];    d=dim(microbiome)[2]
    datasplit<-list();    sid_list=list()
    for (i in 1:length(unique(dat$pop))){ datasplit[[i]]=filter(dat, pop==unique(dat$pop)[i])
      sid_list[[i]]=as.vector(datasplit[[i]]$Sample.ID);      datasplit[[i]]=datasplit[[i]] %>% select(-c(Sample.ID,pop))
      datasplit[[i]]=datasplit[[i]]+1};
    mean_ref_logscale=matrix(0,length(unique(dat$pop)),1);    normalize_logscale=list() 
  for (i in 1:length(unique(dat$pop))){
      mean_ref_logscale[i]=mean(rowMeans(log(as.matrix(datasplit[[i]])),na.rm=T),na.rm=T)
      normalize_logscale[[i]]=rowMeans(log(as.matrix(datasplit[[i]])),na.rm=T)-mean_ref_logscale[i]}
      normalize_gm=cbind(sid_list[[1]],normalize_logscale[[1]])
  #print(normalize_logscale)
  for (i in 2:length(unique(dat$pop))){
    temp=cbind(sid_list[[i]],normalize_logscale[[i]]);      normalize_gm=rbind(normalize_gm,temp)}
    #print(normalize_gm)
    normalize_gm=data.frame(normalize_gm, stringsAsFactors =FALSE); colnames(normalize_gm)=c("Sample.ID","normal")
    #print(normalize_gm$normal)
      return(normalize_gm)}

