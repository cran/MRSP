#include <math.h>
#include <R.h>
#include "mlsp.h"

void maxRow(const double *x, const int *dn, const int *p, double *out){
 int i, j;
 for(i = 0; i < *dn; i++){
  out[i] = x[i];
  for(j = 1; j < *p; j++)
   if(x[i + *dn*j] > out[i])
    out[i] = x[i + *dn*j];
 }
}   

void rowCumsum(const double *x, const int *n, const int *p, double *out){
 double sum;
 int i, j;
 for(i = 0; i < *n; i++){
  sum = 0.0;
  for(j = 0; j < *p; j++){
   sum += x[i + *n*j];
   out[i + *n*j] = sum;
  }
 }
} 

void rowCumprod(const double *x, const int *n, const int *p, double *out){
 double prod;
 int i, j;
 for(i = 0; i < *n; i++){
  prod = 1.0;
  for(j = 0; j < *p; j++){
   prod *= x[i + *n*j];
   out[i + *n*j] = prod;
  }
 }
} 

void yseqlogconstructor(const double *xn, const double *xd, const int *n, const int *p, const int *m, double *out){
 int i, j;
 for(i = 0; i < *n; i++){
  for(j = 0; j < *p; j++){
   if(m[i] < j) out[i + *n*j] = 0.0;
   else if(m[i] == j) out[i + *n*j] = 1.0 / xd[i + *n*j];
   else if(m[i] > j) out[i + *n*j] = xn[i + *n*j] / xd[i + *n*j];
  }
 }
}

void lassoC(const double *u, double *lambda, const double *w, const int *dp, const int *dn, double *norms, double *out){
 int i, j;
 double s;
 if(*dp == 1){
  for(i = 0; i < *dn; i++){
   norms[i] = *lambda * w[i];
   if(u[i] > norms[i]) out[i] = u[i] - norms[i];
   else if(u[i] < -norms[i]) out[i] = u[i] + norms[i];
   else out[i] = 0.0;
  }
 }else{
  *lambda = *lambda * sqrt(*dp);
  for(i = 0; i < *dn; i++){
   for(j = 0; j < *dp; j++){
    norms[i] = norms[i] + u[i + *dn*j] * u[i + *dn*j];
   }
   norms[i] = sqrt(norms[i]);
   if(norms[i] >  *lambda * w[i]){
    s = 1 - *lambda * w[i] / norms[i];
    for(j = 0; j < *dp; j++){
     out[i + *dn*j] = u[i + *dn*j]*s;
    }
   }else{
    for(j = 0; j < *dp; j++){
     out[i + *dn*j] = 0.0;
    }
   }
  }
 }
}  



