
reset_dat_missing <- function(dat,pii=NULL, ref_name=NULL, zeros){
  if(is.null(ref_name)==T){normal_frame=normalizer_gm(dat)}
  if(is.null(ref_name)==F){normal_frame=normalizer_ref(dat,ref_name)}
  if(is.null(pii)==T){pii=0.25}; 
    #if(is.null(ref)==T){ref=nz_ref(dat)}
  ###########ref_name is REFERENCE MICORBE##################################
  #print(dat)
  dat=merge(normal_frame,dat,by="Sample.ID")
  #print(dat)
  datasplit<-list();sid_list=list();pop_list=list(); normal_list=list()
  for (i in 1:length(unique(dat$pop))){datasplit[[i]]=dat[dat$pop==unique(dat$pop)[i],] #datasplit[[i]]=filter(dat, pop==unique(dat$pop)[i])

  sid_list[[i]]= datasplit[[i]]$Sample.ID; pop_list[[i]]=datasplit[[i]]$pop
  normal_list[[i]]= datasplit[[i]]$normal; datasplit[[i]]=datasplit[[i]] %>% select(-c(pop,Sample.ID,normal))}
  popu=unlist(pop_list); sid=unlist(sid_list); normalu=unlist(normal_list)
  sub_dat=datasplit;   sub_dat_adjust=sub_dat
  d=dim(sub_dat[[1]])[2] ;  adj_microbes=matrix(0,length(unique(dat$pop)),d)
  #print(normal_list) 
  for (i in 1:length(unique(dat$pop))){sub_dat_pseudo=sub_dat[[i]]+1 ;    nsub=dim(sub_dat_pseudo)[1]
  d=dim(sub_dat_pseudo)[2];     normalizer=matrix(0,nsub,1)
  normalizer=as.vector(as.numeric(normal_list[[i]]))
  all_ind=which(colSums(zeros)==0);  
  if(is.null(ref_name)==T){ind_adjust=all_ind};  
  if(is.null(ref_name)==F){ind_adjust=setdiff(all_ind, which(names(sub_dat_pseudo)==ref_name))}
  for (j in ind_adjust){
  #print(normalizer)
  microbe=log(as.numeric(sub_dat_pseudo[,j]))-normalizer
  #print("###############\n")
  #print(microbe)
  aa=pop_detect(microbe) ;    cluster_seq=aa[[7]]
  par=c(aa[[1]],aa[[2]],aa[[3]],aa[[4]],aa[[5]]) ;  rightend=par[[1]]+(1.96)*par[[3]]
  leftend=par[[2]]-(1.96)*par[[4]];     pileft=par[[5]];  piright=1-par[[5]];   rml=0;
  if(pileft<pii){rml=1} ;  
  if((rightend<leftend) && (rml==1)){if(length(cluster_seq==0)>0){ 
  adj_microbes[i,j]=1; clus_ind=which(cluster_seq==0) ; sub_dat_adjust[[i]][clus_ind,j]=NA}}}}
  
  
  
  
  number_adj_microbes=list() ;  names_adj_microbes=list()
  
  
  for (i in 1:length(unique(dat$pop))){number_adj_microbes[[i]]=sum(adj_microbes[i,])
  
  
  names_adj_microbes[[i]]=names(datasplit[[i]])[which(adj_microbes[i,]==1)]}
  
  
  
  
  
  dat_adjust=sub_dat_adjust[[1]]
  for (i in 2:length(unique(dat$pop))){dat_adjust=rbind(dat_adjust,sub_dat_adjust[[i]])}
  
  
  
  
  
  
  dat_adjust=data.frame(Sample.ID=sid,pop=popu, normal=normalu,dat_adjust)
  
  
  out=list(dat_adjust, number_adj_microbes, names_adj_microbes,adj_microbes)
  names(out)=c("adjusted data", "number of adjusted microbes", "names of adjusted microbes", "adjusted microbes")
  return(out)}



