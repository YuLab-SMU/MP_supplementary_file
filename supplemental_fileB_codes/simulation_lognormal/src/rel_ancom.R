rel_ancom<-function(OTUdat, Vardat, main.var="country", pr=NULL, pii=NULL,ref_name=NULL){
  source("./src/pop_detect.R");  source("./src/stepA1.R");  source("./src/stepA2.R");  source("./src/geom_ref.R")
  ########## STEP A1 ################################################################
  A1=struc_zero(OTUdat = OTUdat, Vardat = Vardat, main.var = main.var, p=pr)
  #print(A1[[2]])
  ########## STEP A2 ################################################################
  A2=reset_dat_missing(dat=A1[[1]],pii=pii,ref_name=ref_name,zeros=A1[[2]])
  #print(A2[[1]])
  adjusted_data=A2[[1]]; number_A2_adj_microbes=A2[[2]]; names_A2_adj_microbes=A2[[3]]
  ########## STEP A3 ################################################################
  adjusted_data[adjusted_data==0]=1
  #print(adjusted_data)
  ########## STEP B ################################################################
  ########## filter struc zero################################################################
  struc_zero_elim_ind= which(colSums(A1[[2]])>1);   sid=adjusted_data[,1];  
  normalu=as.numeric(adjusted_data[,3]);   popu=adjusted_data[,2];  dat=adjusted_data[,-c(1,2,3)]
  ###Normalize and test individually
  if (length(struc_zero_elim_ind)==0){redat=log(dat)- normalu}else{
  redat=log(dat[,-struc_zero_elim_ind])- normalu}
  gauss_data=data.frame(Sample.ID=sid, pop=popu, redat)
  #############subset only the microbes#######################################
  gauss_microbiome=gauss_data %>% select(-c(Sample.ID,pop))
  ########################analysis ##########################################
  d=dim(gauss_microbiome)[2]; pval=matrix(0,d,1); sig_ind=matrix(0,d,1)
  #print(gauss_microbiome)
  for (j in 1:d){covariate=gauss_data$pop; response=gauss_microbiome[,j]
  regress=data.frame(response,covariate);  regress2=na.omit(regress)
  model=lm(regress2$response~ factor(regress2$covariate));  anv=anova(model);  pval[j]=anv$`Pr(>F)`[[1]]}
  padj=p.adjust(pval,"BH"); ind=which(padj<0.05)
  detected_microbes=names(gauss_microbiome)[ind]
  structural_zeros= data.frame(Microbe=names(adjusted_data[-c(1,2,3)]),t(1-A1[[2]]))
  return(list(detected_microbes=detected_microbes, structural_zeros=structural_zeros))}
####################################################################################
