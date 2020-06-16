module commonpar
  ! for global parameters because integrate (dqags) doesn't know about
  !   other parameters in the function being integrated
  double precision gpar(3)
  integer xx 
end module

!rogram main
! use commonpar
! implicit none
! double precision I2gkxchf
! double precision invchf 
! double precision u,rtem,temi,param(3)
! double precision pmfvec(10),cdfvec(10)
! integer x,ix,xmx
! external I2gkxchf
! param(1)=0.5d0; param(2)=0.1d0; param(3)=4.0d0
! xmx=4
! read *, param(1:3)
! print *, "***checking getcdfpmf for sum_{i=1}^xprev K_i(alp,gam) I2"
! print *, param
! gpar=param
! temi=invchf(x,I2gkxchf,param)
! print *,temi
! call getcdfpmf(xmx,I2gkxchf,param,cdfvec,pmfvec)
! do ix=0,xmx
!   print *, ix,pmfvec(ix+1),cdfvec(ix+1)
! end do
! stop
! end

! front end to link to C (or fortran90)
! extended binomial distribution
! xmx = upper bound for which to get cdf and pmf
! icode=1 for I1 (binomial thinning), icode=2 (family I2), icode=3 (family I3)
! param = (alpha, gamma, xprev)
! output pmfvec(1:(xmx+1)),cdfvec(1:(xmx+1)) allocated in calling routine
subroutine ebinom(xmx,icode,param,cdfvec,pmfvec)
  implicit none
  integer xmx, icode, ix
  double precision I1gkxchf, I2gkxchf, I3gkxchf
  double precision param(3), cdfvec(xmx+1),pmfvec(xmx+1),xprev
  external I1gkxchf, I2gkxchf, I3gkxchf
  select case(icode)
  case(1)
    call getcdfpmf(xmx,I1gkxchf,param,cdfvec,pmfvec)
  case(2)
    call getcdfpmf(xmx,I2gkxchf,param,cdfvec,pmfvec)
  case(3)
    call getcdfpmf(xmx,I3gkxchf,param,cdfvec,pmfvec)
  end select
  return
  end

! xmx = upper bound for which to get cdf and pmf
! fnchf = integrand for characteristic function chf
! param = vector with dimension 3  (alpha,gamma,xprev)
subroutine getcdfpmf(xmx,fnchf,param,cdfvec,pmfvec)
  implicit none
  double precision fnchf, param(3), cdfvec(xmx+1),pmfvec(xmx+1),xprev
  double precision invchf 
  integer xmx, ix
  external fnchf
  xprev=param(3)
  if(xprev>0) then
    do ix=0,xmx
      !print *, "x=", ix
      cdfvec(ix+1)=invchf(ix,fnchf,param)
    end do
    else
    cdfvec(1:(xmx+1))=1    
  endif
  pmfvec(1)=cdfvec(1)
  do ix=1,xmx
    pmfvec(ix+1)=cdfvec(ix+1)-cdfvec(ix)
  end do
  return
  end
  
! invert pgf after it is converted to a characteristic function
!   with integrand fnc and parameter param
! x = non-negative integer
! fnc = function for integrand
! param = parameter vector for fnc
double precision function invchf(x,fnc,param)
  use commonpar
  implicit none
  double precision fnc, param(3)
  integer x
  ! next lines for dqags (integration over bounded interval)
  double precision abserr, work(4000), epsabs, epsrel
  double precision pi,npi, integ1,integ2,integv
  integer x1, neval, ier, limit, lenw, last, iwork(1000)
  external fnc
  epsabs=1.d-6; epsrel=1.d-3
  limit=1000; lenw=limit*4
  pi=3.14159265358979323846d0; npi=-pi
  x1=x+1
  xx=x1; gpar=param

  call dqags(fnc,npi,0.d0,epsabs,epsrel,integ1,abserr,neval,ier,limit,lenw,last,iwork,work)
  call dqags(fnc,0.d0,pi,epsabs,epsrel,integ2,abserr,neval,ier,limit,lenw,last,iwork,work)
  !print *, integ1,integ2
  integv=(integ1+integ2)/(2.d0*pi)
  invchf=0.5d0-integv
  return
  end

! chf for compounding with I3 operator  sum_{j=0}^xprev K_j(x;alpha,gamma)
! xx and param are global variables in commonpar
! u = argument of integrand
double precision function I3gkxchf(u)
  use commonpar
  implicit none
  double precision u,rtem,alp,gamma,xprev
  double complex onei,z,ph,gk,tem
  alp=gpar(1); gamma=gpar(2); xprev=gpar(3)
  onei=cmplx(0.d0,1.d0);
  z=exp(u*onei);
  gk=(1/gamma)*(1+gamma-(1+gamma-gamma*z)**alp) ! chf of I3 operator
  ph=gk**xprev
  tem=ph*exp(-u*xx*onei)/(1.d0-exp(-u*onei))
  rtem=real(tem)
  I3gkxchf=rtem
  !print *,u,z,ph,tem,rtem
  return
  end

! chf for compounding with I2 operator  sum_{j=0}^xprev K_j(x;alpha,gamma)
! xx and param are global variables in commonpar
! u = argument of integrand
double precision function I2gkxchf(u)
  use commonpar
  implicit none
  double precision u,rtem,alp,gam,xprev,a1
  double complex onei,z,ph,gk,tem
  alp=gpar(1); gam=gpar(2); xprev=gpar(3)
  a1=1.d0-alp
  onei=cmplx(0.d0,1.d0);
  z=exp(u*onei);
  gk=(a1+(alp-gam)*z)/(1-alp*gam-a1*gam*z) ! chf of I2 operator 
  ph=gk**xprev
  tem=ph*exp(-u*xx*onei)/(1.d0-exp(-u*onei))
  rtem=real(tem)
  I2gkxchf=rtem
  !print *,u,z,ph,tem,rtem
  return
  end
   
! chf for I1 operator  = binomial thinning, should get binomial distribution
! this is a check for the code
! xx and param are global variables in commonpar
! u = argument of integrand
double precision function I1gkxchf(u)
  use commonpar
  implicit none
  double precision u,rtem,alp,xprev,a1
  double complex onei,z,ph,gk,tem
  alp=gpar(1); xprev=gpar(3)
  a1=1.d0-alp
  onei=cmplx(0.d0,1.d0);
  z=exp(u*onei);
  gk=(a1+alp*z) ! chf of Bernoulli 
  ph=gk**xprev
  tem=ph*exp(-u*xx*onei)/(1.d0-exp(-u*onei))
  rtem=real(tem)
  I1gkxchf=rtem
  !print *,u,z,ph,tem,rtem
  return
  end

