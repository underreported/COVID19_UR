#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* sample main program to call ebinom_ in invchf.f90 */
/* Version for linking with R */

#define NN 100   // an upper bound on length of series
int main(int argc, char *argv[])
{
double params[2]; // alpha, lambda
int i,k,n,x[NN],iprint,icode,xmx; 
double tem;
void ibp(int *n0, int *xvec, double *params, int *iprint, int *icode, double *nllk);
scanf("%d", &n); 
for(i=0;i<n;i++) {scanf("%d", &x[i]);}
iprint=1;
icode=1; 
for(k=0;k<2;k++) {scanf("%lf", &params[k]);} 
ibp(&n,x,params,&iprint,&icode,&tem);
printf("negative log-likelihood = %f\n", tem);
return(0);
}

/* negative log-likelihood at a fixed parameter vector
based on INAR(1) with binomial thinning operator and Poisson innovations*/
/* n0 = legth of series
xvec = count series of length n
params = 2-dimensional vector with (alpha, xprev)  
iprint = 1 for diagnostics/debugging
Returns:
nllk = negative log-likelihood */

void ibp(int *n0, int *xvec, double *params, int *iprint, int *icode, double *nllk)
{ 
double eparams[3];
double *cdfvec,*pmfvec,*innovpmf;
double alpha,lambda,xprob,tem;
int i,j,k,n,y,xmx,gx,gxk,gxi;
double dpois(int x, double lambda);
void ebinom_(int *xmx, int *icode, double *param, double *cdfvec, double *pmfvec);
n= *n0;
alpha=params[0]; lambda=params[1]; 
eparams[0]=alpha;
tem=0;
for(j=0,xmx=0;j<n;j++) {if(xvec[j]>xmx) xmx=xvec[j];}
cdfvec=(double *) malloc((xmx+1) * sizeof(double));
pmfvec=(double *) malloc((xmx+1) * sizeof(double));
innovpmf=(double *) malloc((xmx+1) * sizeof(double));
for(k=0;k<=xmx;k++) {innovpmf[k]=dpois(k,lambda);}
for(i=1;i<n;i++) 
{
eparams[2]=xvec[i-1];
y=xvec[i];
ebinom_(&xmx,icode,eparams,cdfvec,pmfvec);
for(k=0,xprob=0;k<=y;k++) {xprob+=pmfvec[k]*innovpmf[y-k];}
tem-=log(xprob);
    
#ifdef MAIN
if(*iprint >=1) printf("%d %f\n", i,tem);
#endif
}
*nllk=tem;
free(cdfvec); free(pmfvec); free(innovpmf);
}

/* ============================================================ */

double dpois(int x, double lambda)
{ 
double pmf;
if (x<0) return 0.0;
pmf=exp(x*log(lambda)-lambda-lgamma(x+1.0)); 
return(pmf);
}


