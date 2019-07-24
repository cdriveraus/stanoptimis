d<-function(){


library(rstan)

stext<-'
data{
int dim;
int n;
vector[dim] dat[n];
}
parameters{
vector[dim] means;
cholesky_factor_corr[dim] mcholcor;
vector[dim] logscale;
}
model{
dat ~ multi_normal_cholesky(means,diag_pre_multiply(exp(logscale),mcholcor));
means ~ normal(0,1);
mcholcor ~ lkj_corr_cholesky(2);
logscale ~ normal(-1,1);
}
generated quantities{
matrix[dim,dim] mcov = quad_form_diag( mcholcor * mcholcor\', exp(logscale));
}
'

dim=as.integer(20)
n=as.integer(100)
covsqrt <- matrix(rnorm(dim^2),dim,dim)
cholm <- t(chol(covsqrt %*% t(covsqrt)))
covm <- cholm %*% t(cholm)
means <- rnorm(dim)
dat <- t(apply(matrix(rnorm(dim*n),dim,n),2,function(x) cholm %*% x + means))
dim(dat)

 sm<-stan_model(model_code=stext)

 sfit <- sampling(object = sm,data=list(dat=dat,dim=dim,n=n),control=list(max_treedepth=7),chains=3,cores=3)
 e=extract(sfit)
 s=summary(sfit)$summary
 s

 optimis <- stanoptimis(standata = list(dat=dat,dim=dim,n=n),sm = sm,verbose=1,tdf=2,
   finishsamples = 3000,cores=4)
 optimis

apply(optimis$rawposterior,2,mean)
apply(optimis$rawposterior,2,sd)
isdiag(optimis)

for(i in 1:dim){
plot(density(optimis$transformedpars$means[,i]))
points(density(e$means[,i]),type='l',col=2)

for(j in 1:dim){
plot(density(optimis$transformedpars$mcov[,i,j]))
points(density(e$mcov[,i,j]),type='l',col=2)
}}


library(loo)
ps=psis(-optimis$isdiags$lpsamples)
ps
ps$diagnostics
plot(ps)
pareto_k_values(ps)
}
