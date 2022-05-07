#A.	Clustering algorithm

pop_detect <- function(obs){
  TOL = 1e-2;
  MAXITS=1000
  its = 1;
  change = 1;
  obs=obs
  q1=as.numeric(quantile(obs, 0.25))
  q3=as.numeric(quantile(obs, 0.75))
  iqr=q3-q1
  mu1old=q1
  mu2old=q3
  si1old=iqr/2
  si2old=iqr/2
  yold=ifelse((obs<q1+(iqr/2)), 0,1)
  piold=length(which(yold==0))/length(obs)
  while((its <= MAXITS) && (change > TOL) ){
    ynew=ifelse((dnorm(obs, mean=mu1old,sd=si1old)/dnorm(obs, mean=mu2old,sd=si2old))>((1-piold)/piold), 0,1)
    obs_frame=data.frame(obs,ynew)
    obsnew1=filter(obs_frame, obs_frame$ynew==0)[,1]
    obsnew2=filter(obs_frame, obs_frame$ynew==1)[,1]
    pinew=length(which(ynew==0))/length(obs)
    if(pinew>0.999){pinew=0.999}
    if(pinew<0.001){pinew=0.001}
    mu1new=ifelse(length(obsnew1)==0, mean(obsnew2), mean(obsnew1))
    si1new=ifelse(length(obsnew1)==1, 10^(-4),sd(obsnew1))
    si1new=ifelse(length(obsnew1)==0, sd(obsnew2), si1new)
    mu2new=ifelse(length(obsnew2)==0,mean(obsnew1),mean(obsnew2))
    si2new=ifelse(length(obsnew2)==1, 10^(-4), sd(obsnew2))
    si2new=ifelse(length(obsnew2)==0, sd(obsnew1), si2new)
    change = max(abs(mu1new-mu1old),abs(mu2new-mu2old),abs(pinew-piold),abs(si1old-si1new),abs(si2old-si2new)) ;
    mu1old=mu1new
    mu2old=mu2new
    si1old=si1new
    si2old=si2new
    piold=pinew
    its = its + 1;}
  if(its == MAXITS + 1){
    sprintf('max iterations');}
  return(list(mu1new,mu2new,si1new,si2new,pinew, obs, ynew, its));}


