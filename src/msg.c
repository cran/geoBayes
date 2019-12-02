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

void F77_SUB(msgmci2)(int *i, int *a, int *b)
{
  Rprintf("%9d  %8d  %8d\n", *i, *a, *b);
}

void F77_SUB(msgmca2)(void)
{
  Rprintf("Iteration  Acc%% cov    Acc%% z\n");
}

void F77_SUB(msgmce2)(int *a, int *b)
{
  Rprintf("Avg acc %%  %8d  %8d\n", *a, *b);
}

void F77_SUB(msgmcl2)(void)
{
  Rprintf("-----------------------------\n");
}
