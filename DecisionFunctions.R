#requires mixtools


GaussMixBased <-function(dep_data_cols, noise_percents, return_noiseless=FALSE, return_fitted=FALSE, p_0=(781/(781+613)), randseed=NULL)
{   
if(!is.null(randseed)) set.seed(randseed)    
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
current_prop=0.5
for (i in seq_along(IDX)){
  f_1=EMfits[[i]]$lambda[[1]]*dnorm(d4t4[,IDX[[i]]]+NOISES[[1]][,i] ,EMfits[[i]]$mu[[1]], EMfits[[i]]$sigma[[1]]) 
  f_2=EMfits[[i]]$lambda[[2]]*dnorm(d4t4[,IDX[[i]]]+NOISES[[2]][,i] , EMfits[[i]]$mu[[2]], EMfits[[i]]$sigma[[2]])
  if ((sum(f_1>f_2)/length(f_1))>0.5 ) j_max=1 else j_max=2
  j_min=ifelse(j_max==1,2,1)
  if (current_prop < p_0){ j1=j_min; j2=j_max }
                        else{j1=j_max;j2=j_min}
  f_jk=list(f_1,f_2)
  y1[,i]=f_jk[[j1]]
  y2[,i]=f_jk[[j2]]
  J1[i]<-j1
  J2[i]<-j2
  if (i>1)
    current_S=rowSums(log(y1[,1:i]/y2[,1:i]))
  else
    current_S=log(y1[,1:i]/y2[,1:i])
  current_prob= (exp(current_S))/(1+exp(current_S))
  if (sum(is.na(current_prob))!=0)
        { message('NA estimated probability, try again');
          return(NULL)    }
  #print(i)
  #print(y1[,1:i][1:20])
  #print(y2[,1:i][1:20])
  #print(current_S[1:20])
  #print(current_prob[1:20])
  current_prop=sum( (unlist(lapply(current_prob, function(pr) sample(c(0,1),1,prob=c(1-pr,pr)))))==0)/nrow(d4t4)
  print(current_prop)
}

#Y= rowSums(log((y1)/(y2))) >0
#log_p <-function(X)
#{       
#    ifelse(X<=0, 0.0000001,log(X))
#}
S= rowSums(log((y1)/(y2))) 
p_y=exp(S)/(1 +exp(S))
Y= unlist(lapply(p_y, function(pr) sample(c(1,0),1,prob=c(pr,1-pr)))) 
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
    Scl= rowSums(log((y1cl)/(y2cl))) 
    p_ycl=exp(Scl)/(1 +exp(Scl))
    Y_clean= unlist(lapply(p_ycl, function(pr) sample(c(1,0),1,prob=c(pr,1-pr)))) 
   # Y_clean= rowSums(log((y1cl)/(y2cl))) >0
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


ReverseBayes<-function(dep_data_cols, noise_percents,P,quants,  randseed=NULL)
{
if( !is.null(randseed)) set.seed(randseed)    
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

SumOf1s<-function(dep_data_cols, noise_percents,quants,p_k=c(rep(0,4),0.3,0.5,0.7,rep(1,3)),  randseed=NULL)
{
if (!is.null(randseed)) set.seed(randseed)    
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

EM2weighted.density<- function(EM, X, `%noi`=0)
    {
        do.call(cbind,lapply(seq_along(EM$lambda),
                             function(h) 
                                 EM$lambda[[h]]*dnorm(X + rnorm(length(X),0, `%noi`*EM$sigma[[h]]), 
                                                      EM$mu[[h]],
                                                      EM$sigma[[h]])
                             )
                )
    }

expand.grid.list<-function(..., stringsAsFactors=F){
    tdf<-expand.grid(..., stringsAsFactors = stringsAsFactors)
    lapply(split(tdf,seq(nrow(tdf))),as.list)
} 

GaussANDMix<- function(dep_data_cols, dep_data_MMfits, AND.str=0, mg=30, mA=30, AND.dist.from.mu.in.sds=1, noise_percent=0.)
{
    EMfits<- dep_data_MMfits
    X_k<- dep_data_cols
    n_k<- ncol(X_k)
    n_h<-2
    if (n_k<2) stop('at least 2 cols in dep_data_cols are required')
    n_obj<- nrow(X_k)
    gaussians<-expand.grid.list( lapply(1:n_k, function(.dummy) 1:n_h ) )
    gaussians<-lapply(gaussians, unlist)
    w.densities=lapply(1:n_k, function(k) EM2weighted.density(EMfits[[k]], X_k[,k], noise_percent )
                       )
    gaussian.centers<- lapply(gaussians, function(whichMus) 
                              unlist(lapply(seq_along(whichMus),function(h) 
                                    EMfits[[h]]$mu[[ whichMus[[h]] ]] 
                                           )
                                    )
                              )
    gaussian.sigmas<- lapply(gaussians, function(whichMus) 
                              unlist(lapply(seq_along(whichMus),function(h) 
                                    EMfits[[h]]$sigma[[ whichMus[[h]] ]] 
                                           )
                                    )
                              )
    mu.max=apply(do.call(rbind,gaussian.centers),
                 2,max)
    CM<-do.call(rbind,gaussian.centers)
    which.gauss.class1<- which(apply(do.call(cbind,
                                       lapply(1:ncol(CM), function(i) 
                                           CM[,i]==mu.max[[i]]
                                              )
                                        ),
                                     1, sum)==2)
    c1.gauss<- gaussians[[which.gauss.class1]]
    c1.center.dist<- rowSums( (X_k- gaussian.centers[[which.gauss.class1]])^2)^(1/2)
c1.score<-apply(do.call(cbind, 
                  lapply(seq_along(w.densities), function(k) 
                    {
                    w.densities[[k]][,c1.gauss[[k]]]
                    }  )
                 ),
                1,min
               )

AND.score<- apply(X_k +rep(noise_percent,n_k)*gaussian.sigmas[[which.gauss.class1]]  - (mu.max - AND.dist.from.mu.in.sds*gaussian.sigmas[[which.gauss.class1]] ),1,min) 
prg= atan_py(c1.score - median(c1.score),m=mg)
prA=atan_py(AND.score,m=mA)
if( AND.str==0) wA<-0 else
                        wA= ((c1.center.dist)/max(c1.center.dist))^(1/AND.str)
prs=(1-wA)*prg + wA* prA
Y= unlist(lapply(prs, function(pr) sample(c(1,0),1, prob =c(pr,1-pr))  ))
    list(Y=Y,
        w.densities=w.densities,
        fits=EMfits,
        which.density.for.each.variable.is.class1=c1.gauss,
        mu.max=mu.max,
        centers=gaussian.centers,
        sigmas=gaussian.sigmas)
}

