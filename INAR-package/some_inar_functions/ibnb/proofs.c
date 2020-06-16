#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* maybe not all of the include files are needed for MacOS */
/*  Y_k = observed (underreported) X_k is latet bin thinning process
    Y_k<=X_k
gam_k(y[1:k],xk,ik)
 = sum_{x(k-1)>=y(k-1), x(k-1)<=B, i(k-1)=0,1} P(Yk=yk}Xk=xk, Ik=ik)
   * P(Xk=xk|X(k-1)=x(k-1))
   * P(Ik=ik|I(k-1)=i(k-1))
   * gam_{k-1}(y[1:k-1],x(k-1),i(k-1)
   defined for xk>=yk and ik=0 or 1.

gam1 = P(X1=x1) * P(I1=i1) * P(Y1=y1|x1, i1)  for x2>=y1

gam2 = sum_{x1=y1}^B sum_{i1=0}^1 ...
   x2\in {y2,y2+1,...,B} , i2= {0,1}

  B is chosen based on parameters to cover high quantile of
  univariate marginal distribution;
  B is larger than any observed y value,
  for example, B=3*max(yvec)

storing  gam1 is a Bx2 array
         gam2 is a Bx2 array
    ..
         gamT is a Bx2 array
  Actually, just need gamprev and gamnew to avoid storing as Bx2xT 3-dim array
  log (gamT) is the loglik when sum over xT>=yT and iT.
*/

#ifdef MAIN
#define YMX 50   // a large quantile of marginal distribution
#define NN 100   // an upper bound on length of series
main(int argc, char *argv[])
{ double param[5];
  // parameter=(alpha, lambda, omega, q, p01)  
  int i,k,n,y[NN],ymax,B,mult,iprint; 
  double tem;
  double gprev0[YMX],gprev1[YMX],gnew0[YMX],gnew1[YMX];
  void gprob(int k, int *yvec, double *param, int B, 
     double *gprev0, double *gprev1, double *gnew0, double *gnew1);
  void forwprob(int *n0, int *yvec, double *param, int *mult, int *, double *nllk);
  //param[0]=.5; param[1]=1.; param[2]=.5; param[3]=.5; param[4]=.5;
  /* read in data and then call function forwprob for log-likelihood evaluation */
  scanf("%d", &n); // read in sample size = length of series
  for(i=0,ymax=0;i<n;i++) 
  { scanf("%d", &y[i]);  // read in count data values
    if(y[i]>ymax) ymax=y[i];
  }
  mult=3; iprint=1;
  for(k=0;k<5;k++) scanf("%lf", &param[k]); // read in parameters 
  forwprob(&n, y, param, &mult, &iprint, &tem);
  printf("negative log-likelihood = %f\n", tem);
  /* earlier test
  B=mult*ymax;
  for(i=0;i<n;i++) 
  { gprob(i+1,y,param,B,gprev0,gprev1,gnew0,gnew1);
    for(k=0,tem=0;k<=B;k++) 
    { tem=tem+gnew0[k]+gnew1[k];
      gprev0[k]=gnew0[k];
      gprev1[k]=gnew1[k];
    }
    tem=-log(tem);
    printf("%d %f\n", i+1,tem);
  }*/
  return(0);
}
#else
// this include file might be needed for Rprintf
#include <R.h>
#endif

/* negative log-likelihood at a fixed parameter vector
   based on forward probabilities for the latent Markov chain */
/* n0 = legth of series
   yvec = count series of length n
   param = 5-dimensional vector with (alpha, lambda, omega, q, p01)  
   mult = multiple for max(yvec) for upper bound for summation, e.g. 3
   iprint = 1 for diagnostics/debugging
   Returns:
   nllk = negative log-likelihood 
*/
void forwprob(int *n0, int *yvec, double *param, int *mult, int* iprint, double *nllk)
{ 
  double *gprev0,*gprev1,*gnew0,*gnew1,tem;
  int i,k,n,ymax,B; 
  void gprob(int k, int *yvec, double *param, int B, 
     double *gprev0, double *gprev1, double *gnew0, double *gnew1);
  n= *n0;
  for(i=0,ymax=0;i<n;i++) { if(yvec[i]>ymax) ymax=yvec[i]; }
  B=(*mult)*ymax;
  // allocate space for the vectors
  gprev0=(double *) malloc((B+1) * sizeof(double));
  gprev1=(double *) malloc((B+1) * sizeof(double));
  gnew0= (double *) malloc((B+1) * sizeof(double));
  gnew1= (double *) malloc((B+1) * sizeof(double));
  for(i=0;i<n;i++) 
  { gprob(i+1,yvec,param,B,gprev0,gprev1,gnew0,gnew1);
    for(k=0,tem=0;k<=B;k++) 
    { tem=tem+gnew0[k]+gnew1[k];
      gprev0[k]=gnew0[k];
      gprev1[k]=gnew1[k];
    }
    tem=-log(tem);
#ifdef MAIN
    if(*iprint >=1) printf("%d %f\n", i+1,tem);
#else
    if(*iprint >=1) printf("%d %f\n", i+1,tem);
#endif
  }
  *nllk=tem;
  // deallocate space for the vectors
  free(gprev0); free(gprev1); free(gnew0); free(gnew1);
}

// k = integer between 1 and n
// yvec = data vector of counts
// param = 5-dimensional parameter vector
// B = upper bound to use in summations (should be negligble probability
//  for values about B
// gprev0 = prev gamma when i[k]=0
// gprev1 = prev gamma when i[k]=1
// Returns:
// gnew0 = new gamma when i[k]=0, obtained by recursion
// gnew1 = new gamma when i[k]=1
void gprob(int k, int *yvec, double *param, int B, 
   double *gprev0, double *gprev1, double *gnew0, double *gnew1)
{ double alpha, lambda, omega, q, p01,mu, ia2ib[4];
  double Xtransprob(double xa, double xb, double alp, double lambda);
  void Itransprobm(double omega, double p01, double *ia2ib);
  double eprob(int xb, int ib, int y, double q);
  double dpois(int x, double lambda);
  double xprob,iprob,yprob,gpr0,gpr1,gtem,tem0,tem1;
  int y,x,iprev,xprev,yprev;
 
  alpha=param[0]; lambda=param[1]; omega=param[2]; q=param[3]; p01=param[4];
  mu=lambda/(1.-alpha);
  Itransprobm(omega,p01,ia2ib);
  y=yvec[k-1]; 
  for(x=0;x<=B;x++) gnew0[x]=0.;
  for(x=0;x<=B;x++) gnew1[x]=0.;
  // beginning case, start in state 0
  if(k==1) 
  { // gnew0, gnew1
    for(x=y;x<=B;x++) 
    { xprob=dpois(x,mu);   
      // from R code but this might be incorrect
      iprob=ia2ib[0]; // position 1
      yprob=eprob(x,0,y,q);
      gnew0[x]=xprob*iprob*yprob;
      iprob=ia2ib[1]; // position 2
      yprob=eprob(x,1,y,q);
      gnew1[x]=xprob*iprob*yprob;
    }
  }
  else /* recursion */
  { // gnew0, gnew1 
    yprev=yvec[k-2];
    for(x=y;x<=B;x++) 
    { gpr0=0.; gpr1=0.;
      for(xprev=yprev;xprev<=B;xprev++) 
      { xprob=Xtransprob(xprev,x,alpha,lambda);
        for(iprev=0;iprev<=1;iprev++)
        { 
          if(iprev==0) gtem=gprev0[xprev]; else gtem=gprev1[xprev];
          iprob=ia2ib[2*iprev]; // inew=0, position 1 or 3
          yprob=eprob(x,0,y,q); // inew=0
          tem0=xprob*iprob*yprob*gtem;
          // above or below are same, below not needed to prevent roundoff
          //tem0=exp(log(xprob)+log(iprob)+log(yprob)+log(gtem));
          gpr0=gpr0+tem0;
          iprob=ia2ib[2*iprev+1]; // inew=1, position 2 or 4
          yprob=eprob(x,1,y,q);  // inew=1
          tem1=xprob*iprob*yprob*gtem;
          //tem1=exp(log(xprob)+log(iprob)+log(yprob)+log(gtem));
          gpr1=gpr1+tem1;
        }
      }
      gnew0[x]=gpr0; gnew1[x]=gpr1;
    }
  }
}

// ============================================================

double dpois(int x, double lambda)
{ double pmf;
  if (x<0) return 0.0;
  pmf= exp(x*log(lambda) - lambda - lgamma(x+1.0)); 
  return(pmf);
}

double dbinom(int x, int n, double p)
{ double pmf,bincoef;
  if (x < 0 || x > n) return 0;
  bincoef = lgamma(n+1.) - lgamma(x+1.) - lgamma(n-x+1.);
  pmf = exp(bincoef + log(p)*x + (n-x)*log(1-p));
  return pmf;
}

/* functions for the conditional probabilities */

// latent binomial thinning process, Poisson innovation
// P(X(n)=x1 | X(n-1)=x0)
// data=(x(n-1),x(n)) 
// xa = previous x value
// xb = new x value
// alp = alpha parameter 
// lambda = lambda parameter 
// Returns transition probability P(X(new)=xb | X(prev)=xa)
double Xtransprob(double xa, double xb, double alp, double lambda)
{ int kk,xmin;
  double pr;
  double dbinom(int,int,double);
  double dpois(int,double);
  xmin=xa; if(xb<xmin) xmin=xb;
  for(kk=0,pr=0.;kk<=xmin;kk++)
  { pr+= dbinom(kk,xa,alp) * dpois(xb-kk,lambda); }
  return(pr);
  }
/* t.probs=function(data,par) 
{ kk=0:min(data[1],data[2])
  sum(dbinom(kk,data[1],par[1])*dpois(data[2]-kk,par[2]))
}*/

// transition probability matrix (as vector) for latent state
// P(I(n)=i(n) | I(n-1)=i(n-1))
// p=parameter=(alpha, lambda, omega, q, p01)  
// p11=1-p01*(1-omega)/omega
// parameters omega, p01
// Returns  ia2ib = vector of length 4 with p00, p01, p10, p11
void Itransprobm(double omega, double p01, double *ia2ib)
{ double p10;
  ia2ib[0]=1-p01; ia2ib[1]=p01;
  p10=p01*(1-omega)/omega;
  ia2ib[2]=p10; ia2ib[3]=1.-p10;
}
/* t.matrix=function(p) 
{ p10=p[5]*(1-p[3])/p[3]
  matrix(c(1-p[5],p[5],p10,1-p10),byrow=TRUE, ncol=2)
} */

// transition probability for latent state 
// data=(x(n-1),x(n),    i(n-1),i(n)) 
// replace by iprob(ia,ib)=Itransprobm[2*ia+ib]  pos 0,1,2,or3
/*i.probs=function(data,t.matrix)
{ i.aux=NULL
  i.aux[data[3]==0 & data[4]==0]=t.matrix[1,1]
  i.aux[data[3]==0 & data[4]==1]=t.matrix[1,2]
  i.aux[data[3]==1 & data[4]==0]=t.matrix[2,1]
  i.aux[data[3]==1 & data[4]==1]=t.matrix[2,2]
  i.aux
}*/

// P(Y(n)=y|X(n)=x,I(n)=i(n))
// data=(x(n-1),x(n),    i(n-1),i(n), y(n)) 
// p=parameter=(alpha, lambda, omega, q, p01)  
// xb = value of X (not underreported)
// ib = value of state I
// y = value of Y
// q = parameter 
// Returns P(Y(n)=y | X(n)=xb,I(n)=ib))
double eprob(int xb, int ib, int y, double q)
{ double pmf;
  double dbinom(int,int,double);
  if(y>xb) return(0.);
  if(y==xb && ib==0) return(1.);
  if(y<xb && ib==0)  return(0.);
  /* else y<=xb and ib==1 */
  pmf= dbinom(y,xb,q);
  return(pmf);
}
/*e.probs=function(data,p)
{ e.aux=NULL
  if(data[5]> data[2])              e.aux=0
  if(data[5]< data[2] & data[4]==0) e.aux=0
  if(data[5]==data[2] & data[4]==0) e.aux=1
  if(data[5]<=data[2] & data[4]==1) e.aux=dbinom(data[5],data[2],p[4])
  e.aux
}*/
