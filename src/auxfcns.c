#include <R.h>
#include <Rmath.h>
// #include <R_ext/RS.h>

int F77_SUB(isfinite)(double *x)
{
  return R_FINITE(*x);
}


void F77_SUB(rngini)(void)
{
  GetRNGstate();
}


void F77_SUB(rngend)(void)
{
  PutRNGstate();
}


// double F77_SUB(matern)(double *d, double *kappa)
// {
//   if (*d > 0.0)
//     return pow(*d,*kappa) * bessel_k(*d,*kappa,1.) *
//       pow(2.0,1.0-(*kappa)) / gammafn(*kappa);
//   else
//     return 1.0;
// }

double F77_SUB(randnorm)(void)
{
  return rnorm(0,1);
}

double F77_SUB(randunif)(void)
{
  return runif(0,1);
}

double F77_SUB(randt)(double *d)
{
  return rt(*d);
}

double F77_SUB(randgamma) (double *shp)
{
  return rgamma(*shp,1.0);
}

double F77_SUB(randchisq) (double *df)
{
  return rchisq(*df);
}

// d-p-q for t
double F77_SUB(logprobt)(double *q, double *d)
{
  return pt((*q),(*d),1,1);
}

double F77_SUB(logborpt)(double *q, double *d)
{
  return pt((*q),(*d),0,1);
}

double F77_SUB(quantt)(double *p, double *d)
{
  if (*d < 1.0) {
    double p2m1 = 2*expm1(*p) + 1.0; // as p is log(p) p2m1 = 2*p - 1
    if (p2m1 > 0.0)
      return sqrt(qf(p2m1,1.0,*d,1,0));
    else if (p2m1 < 0.0)
      return -sqrt(qf(-p2m1,1.0,*d,1,0));
    else
      return 0.0;
  }
  else
    return qt((*p),(*d),1,1); // This is slow for d < 1
}

double F77_SUB(logpdft)(double *x, double *d)
{
  return dt((*x),(*d),1);
}

// d-p-q for normal
double F77_SUB(logprobnorm)(double *q)
{
  return pnorm((*q),0.0,1.0,1,1);
}

double F77_SUB(logborpnorm)(double *q)
{
  return pnorm((*q),0.0,1.0,0,1);
}

double F77_SUB(quantnorm)(double *p)
{
  return qnorm((*p),0.0,1.0,1,1);
}

double F77_SUB(logpdfnorm)(double *x)
{
  return dnorm((*x),0.0,1.0,1);
}


// d-p-q for logistic
double F77_SUB(logproblogis)(double *q)
{
  return plogis((*q),0.0,1.0,1,1);
}

double F77_SUB(logborplogis)(double *q)
{
  return plogis((*q),0.0,1.0,0,1);
}

double F77_SUB(quantlogis)(double *p)
{
  return qlogis((*p),0.0,1.0,1,1);
}

double F77_SUB(logpdflogis)(double *x)
{
  return dlogis((*x),0.0,1.0,1);
}

// Utility functions
double F77_SUB(flog1p)(double *x)
{
  return log1p((*x));
}

double F77_SUB(fexpm1)(double *x)
{
  return expm1((*x));
}

double F77_SUB(flog1pexp)(double *x)
{
  return -plogis((*x),0.0,1.0,0,1);
}

double F77_SUB(flog1mexp)(double *x)
{
  //return (*x) - qlogis((*x),0.0,1.0,1,1);
  double negln2 = -0.69314718055994530941723212145817656807550013436025525412;
  if ((*x) >= negln2) {
    return log(-expm1(*x));
  } else {
    return log1p(-exp(*x));
  }
}

double F77_SUB(flogexpm1)(double *x)
{
  return qlogis(-(*x),0.0,1.0,0,1);
}

double F77_SUB(fbesselkratio)(double *x, double *ktop, double *kbot)
{
  if (fabs(*ktop) == fabs(*kbot)) return 1.0;
  // The scaling cancels in the ratio.
  return bessel_k(*x,*ktop,2.)/bessel_k(*x,*kbot,2.);
}

double F77_SUB(flogbesselkdk)(double *x, double *kappa)
{
  // Numerical derivative of log(BesselK) w.r.t order for positive order
  double eps = sqrt(DBL_EPSILON);
  // The scaling cancels in the difference.
  if (*kappa < eps)
    return (log(bessel_k(*x,*kappa+eps,2.)) - log(bessel_k(*x,*kappa,2.)))/eps;
  return .5*(log(bessel_k(*x,*kappa+eps,2.)) - log(bessel_k(*x,*kappa-eps,2.)))/eps;
}

double F77_SUB(fbesselk)(double *x, double *nu)
{
  return bessel_k(*x, *nu, 1.);
}

double F77_SUB(fbesselkexp)(double *x, double *nu)
{
  return bessel_k(*x, *nu, 2.);
}

double F77_SUB(fgamma)(double *x)
{
  return gammafn(*x);
}

double F77_SUB(fdigamma)(double *x)
{
  return digamma(*x);
}

double F77_SUB(ftrigamma)(double *x)
{
  return trigamma(*x);
}
