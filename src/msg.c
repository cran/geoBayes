#include <R.h>
#include <R_ext/Print.h>

void F77_SUB(msgmci)(int *i, int *a)
{
  Rprintf("%9d  %8d\n", *i, *a);
}

void F77_SUB(msgmca)(void)
{
  Rprintf("Iteration  Acc%% cov\n");
}

void F77_SUB(msgmce)(int *a)
{
  Rprintf("Avg acc %%  %8d\n", *a);
}

void F77_SUB(msgmcl)(void)
{
  Rprintf("-------------------\n");
}
