#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

#define M 2
#define N 3

void multAx(int m,int n,double A[m][n], double X[n],double B[m]){
  bzero(B,m*sizeof(double));
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      B[i]=B[i]+A[i][j]*X[j];
    }
  }
}

void matrixMult(int m, int n, int q, double A[m][n], double B[n][q], double result[m][q]) {
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < q; j++) {
      result[i][j] = 0;
      for(int k = 0; k < n; k++) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

void transposeMatrix(int m, int n, double A[m][n], double result[n][m]) {
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      result[j][i] = A[i][j];
    }
  }
}

void matPrint(int m, int n, double A[m][n]) {
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      printf("%g ", A[i][j]);
    }
    printf("\n");
  }
}

void vecprint(int m,double X[m]){
  matPrint(m,1,(double (*)[1])X);
}

void plufact(int n,double A[n][n],double *P[n]){
  for(int i=0;i<n;i++){
    P[i]=&A[i][0];
  }

  for(int i=0;i<n-1;i++){
    for(int j=i+1;j<n;j++){
      if(fabs(P[i][i])<fabs(P[j][i])){
        double *t=P[i]; P[i]=P[j]; P[j]=t;
      }
    }
    for(int j=i+1;j<n;j++){
      double alpha=P[j][i]/P[i][i];
      for(int k=i+1;k<n;k++){
        P[j][k]-=alpha*P[i][k];
      }
      P[j][i]=alpha;
    }
  }
}

void plusolve(int n,double LU[n][n],double *P[n], double x[n],double b[n]){
  double y[n];
  for(int i=0;i<n;i++){  // Ly=b
    y[i]=b[(P[i]-&LU[0][0])/n];
    for(int j=0;j<i;j++){
      y[i]-=P[i][j]*y[j];
    }
  }

  for(int i=n-1;i>=0;i--){  // Ux=y
    x[i]=y[i];
    for(int j=i+1;j<n;j++){
      x[i]-=P[i][j]*x[j];
    }
    x[i]/=P[i][i];
  }
}

void matinv(int n,double A[n][n],double B[n][n]){
  double LU[n][n],*P[n];
  memcpy(LU,A,sizeof(double)*n*n);
  plufact(n,LU,P);
  double e[n];
  bzero(e,sizeof(double)*n);

  for(int j=0;j<n;j++){
    e[j]=1;
    double v[n];
    plusolve(n,LU,P,v,e);
    for(int i=0;i<n;i++) B[i][j]=v[i];
    e[j]=0;
  }
}

void undsolve(int m, int n, double A[m][n], double x[n], double b[m]) {
  double A_T[n][m];
  transposeMatrix(m, n, A, A_T);
  double A_mult_A_T[m][m];
  matrixMult(m, n, m, A, A_T, A_mult_A_T);
  double A_mult_A_T_inv[m][m];
  matinv(m, A_mult_A_T, A_mult_A_T_inv);
  double A_T_mult_inv[n][m];
  matrixMult(n, m, m, A_T, A_mult_A_T_inv, A_T_mult_inv);
  multAx(n, m, A_T_mult_inv, b, x);
}

int main() {
  double A[M][N] = {
    {1, 2, 3},
    {4, 5, 6},
  };
  double x[N];
  double b[M] = {6, 9};
  double b_check[M];
  undsolve(M, N, A, x, b);
  printf("x=:\n");
  vecprint(N, x);

  printf("Checking computation of x where b check =:\n");
  multAx(M, N, A, x, b_check);
  vecprint(M, b_check);
  return 0;
}
