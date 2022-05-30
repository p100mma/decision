#requires mixtools


GaussMixBased <-function(dep_data_cols, noise_percents, return_noiseless=FALSE, return_fitted=FALSE, randseed=123)
{   
set.seed(randseed)    
d4t4<-dep_data_cols
IDX=c(1:ncol(dep_data_cols))
y1=matrix(nrow=nrow(d4t4),ncol=length(IDX)); y2=matrix(nrow=nrow(d4t4),ncol=length(IDX)); #model densities 
EMfits<- lapply(IDX, function(x) normalmixEM(d4t4[,x]))
for (i in seq_along(IDX)){
  y1[,i]= EMfits[[i]]$lambda[[1]]*dnorm(d4t4[,IDX[[i]]], EMfits[[i]]$mu[[1]], EMfits[[i]]$sigma[[1]])
  y2[,i]= EMfits[[i]]$lambda[[2]]*dnorm(d4t4[,IDX[[i]]], EMfits[[i]]$mu[[2]], EMfits[[i]]$sigma[[2]])
}
if (length(noise_percents)==1) noise_percents=rep(noise_percents,length(IDX))
ranges=lapply(IDX, function(x) range(d4t4[,x]))
noise_amounts= noise_percents*unlist(lapply(ranges, function(x) x[[2]] - x[[1]]))
noises=lapply(1:nrow(d4t4),function(x) unlist(lapply(noise_amounts, function(x) rnorm(1,0,x))))
noises=Reduce(rbind,noises)
Y= rowSums(log((y1+abs(noises))/(y2+abs(noises)))) >0
result_list=list(Y=Y)
if (return_noiseless)
{
    Y_clean=rowSums(log((y1)/(y2))) >0
    result_list$Y_clean=Y_clean
}
if (return_fitted)
{
    result_list$y1=y1
    result_list$y2=y2
    result_list$EMfits=EMfits
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
