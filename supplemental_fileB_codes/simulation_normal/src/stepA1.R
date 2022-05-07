####STEP A1
struc_zero <- function(OTUdat = otu, Vardat = meta, main.var = "population", p=NULL){
  if(is.null(p)==T){p=0.05}
main_var=main.var ;  meta=Vardat
ind=c(which(names(meta)=="Sample.ID"),which(names(meta)==main_var))
red_map=meta[,ind]
merged=merge(red_map, OTUdat, by="Sample.ID")
names(merged)[2]="pop"
#  merged_implement=merged[,-which(names(merged)=="Sample.ID")]
full_dat=merged;  datasplit<-list(); sid_list=list()
  for (i in 1:length(unique(full_dat$pop))){datasplit[[i]]=filter(full_dat, pop==unique(full_dat$pop)[i])
    sid_list[[i]]= datasplit[[i]]$Sample.ID; datasplit[[i]]=datasplit[[i]] %>% select(-Sample.ID)}
#print(datasplit)
#cat("###############\n")
#print(sid_list)
sid=unlist(sid_list)
d=dim(datasplit[[1]])[2]
zeros=matrix(0,length(unique(full_dat$pop)),(d-1))
#print(zeros)
  for (i in 1:length(unique(full_dat$pop))){for (j in 2:d){
      microbe1=datasplit[[i]][,j];  microbe=na.omit(microbe1)
      x=sum(microbe!=0);    n=length(microbe) ;  r=x/n;      
      if(r<=p){datasplit[[i]][,j]=NA};  if(r<=p){zeros[i,(j-1)]=1}}}
dat_adjust_struc=datasplit[[1]]
  for (i in 2:length(unique(full_dat$pop))){dat_adjust_struc=rbind(dat_adjust_struc,datasplit[[i]])}
dat_adjust=data.frame(Sample.ID=sid,dat_adjust_struc)
rownames(zeros)=unique(full_dat$pop)
return(list(dat_adjust,zeros))}
