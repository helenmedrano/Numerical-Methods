#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <strings.h>

extern void matprint(int m,int n,double A[m][n]);
extern void vecprint(int m,double X[m]);
extern void multAx(int m,int n,double A[m][n],
    double X[n],double B[m]);
extern void multATx(int m,int n,double A[m][n],
    double X[m],double B[n]);
extern void cmatprint(int m,int n,complex A[m][n]);
extern void cvecprint(int m,complex X[m]);
extern void cmultAx(int m,int n,complex A[m][n],
    complex X[n],complex B[m]);
extern void LUfact(int n,double A[n][n]);
extern void backsub(int n,double U[n][n],
    double x[n],double b[n]);
extern void LUsolve(int n,double LU[n][n],
    double x[n],double b[n]);
extern void PLUsolve(int n,double LU[n][n],double *P[n],
    double x[n],double b[n]);
extern void cPLUsolve(int n,complex LU[n][n],complex *P[n],
    complex x[n],complex b[n]);
extern void PLUfact(int n,double A[n][n],double *P[n]);
extern void cPLUfact(int n,complex A[n][n],complex *P[n]);
extern double vecnorm1(int n,double x[n]);
extern double matnorm1(int n,double A[n][n]);
extern double vecnorm2(int n,double x[n]);
extern double cvecnorm2(int n,complex x[n]);
extern double matnorm2(int n,double A[n][n]);
extern complex cdotprod(int n,complex x[n],complex y[n]);
extern void cmultMatByConst(int m, int n, complex constant, complex A[n][n]);
extern void cmatrixMult(int m, int n, int q, complex A[m][n], complex B[n][q], complex result[m][q]);
extern void cidentityMatrix(int n, complex A[n][n]);
extern void cmatrixSub(int n, complex A[n][n], complex B[n][n], complex result[n][n]);
