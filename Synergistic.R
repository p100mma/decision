
ConjQuantile <-function(dep_data_cols, noise_percents, offsets, quant)
{
d4t4<-dep_data_cols
IDX=c(1:ncol(dep_data_cols))
if (length(noise_percents)==1) noise_percents=rep(noise_percents,length(IDX))
ranges=lapply(IDX, function(x) range(d4t4[,x]))
noise_amounts= noise_percents*unlist(lapply(ranges, function(x) x[[2]] - x[[1]]))
noises=lapply(1:nrow(d4t4),function(x) unlist(lapply(noise_amounts, function(x) rnorm(1,0,x))))
noises=Reduce(rbind,noises)
q_j=apply(d4t4[,IDX],2,quantile,quant)
`[x_j>q_j]`=t(apply(d4t4[,IDX]+noises,1, function(row) row> q_j) )
Y=  (rowSums(`[x_j>q_j]`)==1)
return(Y)
}

atan_py<- function(x,m)
{
((1/pi)*atan(m*x)) + (1/2)
}

XorFuzzy <-function(dep_data_cols, offsets,  py_fun, py_fun_args)
{
if (ncol(dep_data_cols)!=2) stop("please use 2 variables only.")
x10=x1-offset1 
x2o=x2-offset2
conj1=pmin(-x1o,x2o); conj2=pmin(x1o,-x2o)
xXOR=pmax(conj1,conj2)
py_fun_args$x= xXOR
py=do.call(py_fun, py_fun_args)
Y= unlist( lapply(py, function(pr) sample(c(0,1),1,prob=c(pr,1-pr)) ) ) 
return(Y)
}
