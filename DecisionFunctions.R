#requires mixtools


GaussMixBased <-function(dep_data_cols, noise_percents, return_noiseless=FALSE, return_fitted=FALSE, p_0=(780/(780+620)), randseed=123)
{   
set.seed(randseed)    
d4t4<-dep_data_cols
IDX=c(1:ncol(dep_data_cols))
EMfits<- lapply(IDX, function(x) normalmixEM(d4t4[,x]))
if (length(noise_percents)==1) noise_percents=rep(noise_percents,length(IDX))
#ranges=lapply(IDX, function(x) range(d4t4[,x]))
#noise_amounts= noise_percents*unlist(lapply(ranges, function(x) x[[2]] - x[[1]]))
noise_amounts1= noise_percents*unlist(lapply(EMfits, function(m) m$sigma[[1]]))
noise_amounts2= noise_percents*unlist(lapply(EMfits, function(m) m$sigma[[2]]))
noises1=lapply(1:nrow(d4t4),function(x) unlist(lapply(noise_amounts1, function(x) rnorm(1,0,x))))
noises1=Reduce(rbind,noises1)
noises2=lapply(1:nrow(d4t4),function(x) unlist(lapply(noise_amounts2, function(x) rnorm(1,0,x))))
noises2=Reduce(rbind,noises2)
NOISES=list(noises1,noises2)
y1=matrix(nrow=nrow(d4t4),ncol=length(IDX)); y2=matrix(nrow=nrow(d4t4),ncol=length(IDX)); #model densities 
J1<-vector(length=length(IDX))
J2<-vector(length=length(IDX))
for (i in seq_along(IDX)){
  j1<-sample(c(0,1),1, prob=c(1-p_0,p_0))
  j2<-(!j1)*1.0
  j1<-j1+1
  j2<-j2+1
  J1[[i]]<-j1
  J2[[i]]<-j2
  y1[,i]= EMfits[[i]]$lambda[[j1]]*dnorm(d4t4[,IDX[[i]]]+NOISES[[j1]][,i] ,EMfits[[i]]$mu[[j1]], EMfits[[i]]$sigma[[j1]]) 
  y2[,i]= EMfits[[i]]$lambda[[j2]]*dnorm(d4t4[,IDX[[i]]]+NOISES[[j2]][,i] , EMfits[[i]]$mu[[j2]], EMfits[[i]]$sigma[[j2]])
}

Y= rowSums(log((y1)/(y2))) >0
#log_p <-function(X)
#{       
#    ifelse(X<=0, 0.0000001,log(X))
#}
#S= rowSums(log((y1)/(y2))) 
#p_y=exp(S)/(1 +exp(S))
#Y= unlist(lapply(p_y, function(pr) sample(c(1,0),1,prob=c(pr,1-pr)))) 
result_list=list(Y=Y)
if (return_noiseless)
{ 
    y1cl=matrix(nrow=nrow(d4t4),ncol=length(IDX)); y2cl=matrix(nrow=nrow(d4t4),ncol=length(IDX)); #model densities 
    for (i in seq_along(IDX)){
        j1<-J1[[i]]
        j2<-J2[[i]]
        y1cl[,i]= EMfits[[i]]$lambda[[j1]]*dnorm(d4t4[,IDX[[i]]], EMfits[[i]]$mu[[j1]], EMfits[[i]]$sigma[[j1]])
        y2cl[,i]= EMfits[[i]]$lambda[[j2]]*dnorm(d4t4[,IDX[[i]]], EMfits[[i]]$mu[[j2]], EMfits[[i]]$sigma[[j2]])}
   # print(sum(is.na(y1cl))) 
   # print(sum(is.na(y1)) )
   # print(sum(is.na(y2cl))) 
   # print(sum(is.na(y2)) )
   # Scl= rowSums(log((y1cl)/(y2cl))) 
   # p_ycl=exp(Scl)/(1 +exp(Scl))
   # Y_clean= unlist(lapply(p_ycl, function(pr) sample(c(1,0),1,prob=c(pr,1-pr)))) 
    Y_clean= rowSums(log((y1cl)/(y2cl))) >0
    result_list$Y_clean=Y_clean
}
if (return_fitted)
{
    result_list$y1=y1
    result_list$y2=y2
    result_list$EMfits=EMfits
    result_list$J1=J1
    result_list$J2=J2
    if (return_noiseless){
        result_list$y1clean=y1cl
        result_list$y2clean=y2cl
        }
}
return(result_list)
}


ReverseBayes<-function(dep_data_cols, noise_percents,P,quants,  randseed=123)
{
set.seed(randseed)    
d4t4<-dep_data_cols
IDX=c(1:ncol(dep_data_cols))
if (length(noise_percents)==1) noise_percents=rep(noise_percents,length(IDX))
ranges=lapply(IDX, function(x) range(d4t4[,x]))
noise_amounts= noise_percents*unlist(lapply(ranges, function(x) x[[2]] - x[[1]]))
noises=lapply(1:nrow(d4t4),function(x) unlist(lapply(noise_amounts, function(x) rnorm(1,0,x))))
noises=Reduce(rbind,noises)
q_j=apply(d4t4[,IDX],2,quantile,quants)
`[x_j>q_j]`=t(apply(d4t4[,IDX]+noises,1, function(row) row> q_j) )
p_j=P*`[x_j>q_j]`
p_j=p_j + ((1-P)*(!`[x_j>q_j]`))
S_j= log (p_j/(1-p_j))
S=rowSums(S_j)
p_y=exp(S)/(1+exp(S))
Y= unlist(lapply(p_y, function(pr) sample(c(1,0),1,prob=c(pr,1-pr))) )
return(Y)
}

SumOf1s<-function(dep_data_cols, noise_percents,quants,p_k=c(rep(0,4),0.3,0.5,0.7,rep(1,3)),  randseed=123)
{
set.seed(randseed)    
d4t4<-dep_data_cols
IDX=c(1:ncol(dep_data_cols))
if (length(noise_percents)==1) noise_percents=rep(noise_percents,length(IDX))
ranges=lapply(IDX, function(x) range(d4t4[,x]))
noise_amounts= noise_percents*unlist(lapply(ranges, function(x) x[[2]] - x[[1]]))
noises=lapply(1:nrow(d4t4),function(x) unlist(lapply(noise_amounts, function(x) rnorm(1,0,x))))
noises=Reduce(rbind,noises)
q_j=apply(d4t4[,IDX],2,quantile,quants)
`[x_j>q_j]`=t(apply(d4t4[,IDX]+noises,1, function(row) row> q_j) )
Sc= rowSums(`[x_j>q_j]`)
p_y= unlist(lapply(Sc, function(s)  if (s>0){  p_k[[s]]} else 0) )
Y= unlist(lapply(p_y, function(pr) sample(c(1,0),1,prob=c(pr,1-pr))))
return(Y)
} 
