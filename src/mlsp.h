#include <Rmath.h>
#include <Rdefines.h>

void maxRow(const double *x, const int *dn, const int *p, double *out);

void rowCumsum(const double *x, const int *n, const int *p, double *out);

void rowCumprod(const double *x, const int *n, const int *p, double *out);

void yseqlogconstructor(const double *xn, const double *xd, const int *n, const int *p, const int *m, double *out);

void lassoC(const double *u, double *lambda, const double *w, const int *dp, const int *dn, double *norms, double *out);
