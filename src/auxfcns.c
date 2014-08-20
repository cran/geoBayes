#include <R.h>
#include <Rmath.h>
// #include <R_ext/RS.h>

void F77_SUB(rngini)(void)
{
  GetRNGstate();
}


void F77_SUB(rngend)(void)
{
  PutRNGstate();
}


double F77_SUB(matern)(double *d, double *kappa)
{
  if (*d > 0.0)
    return pow(*d,*kappa) * bessel_k(*d,*kappa,1) *
      pow(2.0,1.0-(*kappa)) / gammafn(*kappa);
  else
    return 1.0;
}

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
double F77_SUB(logprobt1)(double *q, double *d)
{
  return pt((*q),(*d),1,1);
}

double F77_SUB(logprobt0)(double *q, double *d)
{
  return pt((*q),(*d),0,1);
}

double F77_SUB(quantt)(double *p, double *d)
{
  return qt((*p),(*d),1,1);
}

double F77_SUB(logpdft)(double *x, double *d)
{
  return dt((*x),(*d),1);
}

// d-p-q for normal 
double F77_SUB(logprobnorm1)(double *q)
{
  return pnorm((*q),0.0,1.0,1,1);
}

double F77_SUB(logprobnorm0)(double *q)
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
double F77_SUB(logproblogis1)(double *q)
{
  return plogis((*q),0.0,0.6458,1,1);
}

double F77_SUB(logproblogis0)(double *q)
{
  return plogis((*q),0.0,0.6458,0,1);
}

double F77_SUB(quantlogis)(double *p)
{
  return qlogis((*p),0.0,0.6458,1,1);
}

double F77_SUB(logpdflogis)(double *x)
{
  return dlogis((*x),0.0,0.6458,1);
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

double F77_SUB(flog1mexp)(double *x)
{
  return (*x) - qlogis((*x),0.0,1.0,1,1);
}

double F77_SUB(flogexpm1)(double *x)
{
  return qlogis(-(*x),0.0,1.0,0,1);
}

