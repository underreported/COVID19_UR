# R interface to the C code for evaluating the negative log-likelihood
# of the INAR(1) model with binomial thinning operator and Poisson innovations

# param=(alpha, lambda) 
# xdata = data vector
# iprint >=1 for printing of intermediate calculations
# Returns negative log-likelihood based on forward probabilities

setwd("/Users/Amanda/Dropbox/PhD/papers_by_folders/Paper6-rPackage/some_inar_functions/inAR1")
ibp.nllk=function(param,xdata,iprint=1)
{ 
if(any(param<0)) return(1.e10);
if(param[1]>=1  ) return(1.e10)
if(!is.loaded("ibp")) dyn.load("./ibp.so")
n=length(xdata)
xdata=as.integer(xdata)
if(min(xdata)<0) stop("negative counts")
nllk=0
icode=1
out=.C("ibp", as.integer(n), as.integer(xdata), 
as.double(param), as.integer(iprint),as.integer(icode),
nllk=as.double(nllk))
out$nllk
}

#======================================================================

xdata=c(1,3,0,4,6,7,8,3,3,4)
params=c(0.5,9)
res=ibp.nllk(param=params,xdata=xdata,iprint=2); res
nlm(ibp.nllk,p=params,xdata=data,iprint=1,hessian=TRUE)
