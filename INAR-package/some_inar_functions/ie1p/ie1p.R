# R interface to the C code for evaluating the negative log-likelihood
# of the INAR(1) model with generalized thinning operator (family 1) and Poisson innovations

# param=(alpha,gamma,lambda) 
# xdata = data vector
# iprint >=1 for printing of intermediate calculations
# Returns negative log-likelihood based on forward probabilities
setwd("/Users/Amanda/Dropbox/Thesis/Paper6-rPackage/ie1p")
ie1p.nllk=function(param,xdata,iprint=1)
{ 
if(any(param<=0)) return(1.e10)
if(param[1]>=1  ) return(1.e10)
if(!is.loaded("ie1p")) dyn.load("./ie1p.so")
n=length(xdata)
xdata=as.integer(xdata)
if(min(xdata)<0) stop("negative counts")
icode=2
nllk=0
out=.C("ie1p", as.integer(n), as.integer(xdata), 
as.double(param), as.integer(iprint),
as.integer(icode),nllk=as.double(nllk))
out$nllk
}

#======================================================================

xdata=c(1,3,0,4,6,7,8,3,3,4)
params=c(0.5,0.6,3)
res=ie1p.nllk(param=params,xdata=xdata,iprint=2); res
