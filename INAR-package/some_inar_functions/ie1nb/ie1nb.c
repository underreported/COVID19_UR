#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* sample main program to call ebinom_ in I2invchf.f90 */

#define NN 100   // an upper bound on length of series
int main(int argc, char *argv[])
{
double params[4]; // (alpha, p, theta, gamma)
int i,k,n,x[NN],iprint,icode,xmx; 
double tem;
void ie1nb(int *n0, int *xvec, double *params, int *iprint, int *icode, double *nllk);
scanf("%d", &n); 
for(i=0;i<n;i++) {scanf("%d", &x[i]);}
iprint=1;
icode=2;
for(k=0;k<4;k++) {scanf("%lf", &params[k]);} 
ie1nb(&n, x, params, &iprint, &icode, &tem);
printf("negative log-likelihood = %f\n", tem);
return(0);
}

/* negative log-likelihood at a fixed parameter vector
based on INAR(1) with expectation thinning operator an NB innovations*/
/* n0 = legth of series
xvec = count series of length n
eparams = 3-dimensional vector with (alpha, gamma, xprev)  
iprint = 1 for diagnostics/debugging
Returns:
nllk = negative log-likelihood */

void ie1nb(int *n0, int *xvec, double *params, int *iprint, int *icode, double *nllk)
{ 
double eparams[3];
double *cdfvec,*pmfvec,*innovpmf;
double alpha,p,theta,gamma,xprob,tem;
int i,j,k,n,y,xmx,gx,gxk,gxi;
double dnbinom(int x, int theta, double p);
void ebinom_(int *xmx, int *icode, double *param, double *cdfvec, double *pmfvec);

n= *n0;

alpha=params[0]; p=params[1]; theta=params[2]; gamma=params[3];
eparams[0]=alpha; eparams[1]=gamma;
tem=0;
for(j=0,xmx=0;j<n;j++) {if(xvec[j]>xmx) xmx=xvec[j];}
cdfvec=(double *) malloc((xmx+1) * sizeof(double));
pmfvec=(double *) malloc((xmx+1) * sizeof(double));
innovpmf=(double *) malloc((xmx+1) * sizeof(double));
for(k=0;k<=xmx;k++) {innovpmf[k]=dnbinom(k,theta,p);}
for(i=1;i<n;i++) 
{
eparams[2]=xvec[i-1];
y=xvec[i];
ebinom_(&xmx,icode,eparams,cdfvec,pmfvec);
for(k=0,xprob=0;k<=y;k++) {xprob+=pmfvec[k]*innovpmf[y-k];}
tem=tem-log(xprob);
    
#ifdef MAIN
    if(*iprint >=1) printf("%d %f\n", i,tem);
#endif
}
free(cdfvec); free(pmfvec); free(innovpmf);
*nllk=tem;
}

/* ============================================================ */

double dpois(int x, double lambda)
{ 
double pmf;
if (x<0) return 0.0;
pmf=exp(x*log(lambda)-lambda-lgamma(x+1.0)); 
return(pmf);
}

double dnbinom(int x, int r, double p)
{ 
double pmf,nbincoef;
if (x<0) return 0;
nbincoef=lgamma(r+x)-lgamma(r)-lgamma(x+1.);
pmf=exp(nbincoef) * pow(p,r) * pow((1-p),x);
return(pmf);
}


