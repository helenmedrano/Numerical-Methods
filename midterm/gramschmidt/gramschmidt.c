#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "matrixlib.h"
#include "tictoc.h"

void identityMatrix(int n, double A[n][n]) {
  for(int i = 0; i < n; i++) {
    for(int j  = 0; j < n; j++) {
      if(i == j)
        A[i][j] = 1;
      else
        A[i][j] = 0;
    }
  }
}

void PLUfact(int n, double A[n][n], double P[n][n]) {
  identityMatrix(n, P);
  double Umax, Uii, tmp;
  int q, row, i, r, j;

  for (i=0; i<n; i++) {
    Umax = 0;
    for (r=i; r<n; r++) {
      Uii=A[r][i];
      q = 0;
      while (q<i) {
        Uii -= A[r][q]*A[q][r];
        q++;
      }
      if (fabs(Uii)>Umax) {
        Umax = fabs(Uii);
        row = r;
      }
    }
    if (i!=row) {
      for (q=0; q<n; q++) {
        tmp = P[i][q];
        P[i][q]=P[row][q];
        P[row][q]=tmp;
        tmp = A[i][q];
        A[i][q]=A[row][q];
        A[row][q]=tmp;
      }
    }
    
    j = i;
    while (j<n) { //Determine U across row i
      q = 0;
      while (q<i) {
        A[i][j] -= A[i][q]*A[q][j];
        q++;
      }
      j++;
    }
    j = i+1;
    while (j<n) {
      q = 0;
      while (q<i) {
        A[j][i] -= A[j][q]*A[q][i];
        q++;
      }
      A[j][i] = A[j][i]/A[i][i];
      j++;
    }
  }
  /*
   *   printf("Matrix LU:\n");
   *     matPrint(n, n, A);
   *       printf("Matrix P:\n");
   *         matPrint(n, n, P);*/
}

void matinv(int n,double A[n][n],double B[n][n]){
  double LU[n][n], P[n][n];
  memcpy(LU,A,sizeof(double)*n*n);
  PLUfact(n, LU, P);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      if (P[i][j] == 1)
        B[i][j] = 1.0;
      else
        B[i][j] = 0.0;
      for (int k = 0; k < i; k++)
        B[i][j] -= LU[i][k] * B[k][j];
    }
    for (int i = n - 1; i >= 0; i--) {
      for (int k = i + 1; k < n; k++) {
        B[i][j] -= LU[i][k] * B[k][j];
      }
      B[i][j] = B[i][j] / LU[i][i];
    }
  }
}

double dotprod(int n,double x[n],double y[n]){
    double s=0.0;
    for(int i=0;i<n;i++){
        s+=x[i]*y[i];
    }
    return s;    
}
void gramschmidt(int m, int n,
    double A[m][n],double QT[n][m]){
    double vt[m];
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++) vt[i]=A[i][j];
        for(int k=0;k<j;k++){
            double dp=dotprod(m,vt,QT[k]);
            for(int i=0;i<m;i++){
                vt[i]-=dp*QT[k][i];
            }
        }
        double vtnorm=vecnorm2(m,vt);
        for(int i=0;i<m;i++) QT[j][i]=vt[i]/vtnorm;
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

#define N 5
double A[N][N],QT[N][N],R[N][N];
int main(){
  double b[N];
  for(int i = 0; i < N; i++) {
    b[i] = i;
  }
  double c[N];
//    setrlimit(RLIMIT_STACK,
//        &(const struct rlimit)
//        {RLIM_INFINITY,RLIM_INFINITY});
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
//            A[i][j]=2.0*random()/RAND_MAX-1.0;
            A[i][j]=1.0/(i+j+1);
        }
    }
    printf("N=%d\n",N);
    printf("A=\n"); matprint(N,N,A);
    gramschmidt(N,N,A,QT);
    bzero(R,sizeof(double)*N*N);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                R[i][j]+=QT[i][k]*A[k][j];
            }
        }
    }
    printf("QT=\n"); matprint(N,N,QT);
    printf("R=\n"); matprint(N,N,R);

    printf("b =:\n");
    vecprint(N, b);

    // c = QTb
    matrixMult(N, N, 1, QT, b, c);
    printf("c =:\n");
    vecprint(N, c);

    // Solve Rx = c
    // x = R_invc
    double x[N];
    double R_inv[N][N];
    matinv(N, R, R_inv);
    matrixMult(N, N, 1, R_inv, c, x);
    printf("x solutions =:\n");
    vecprint(N, x);
    


    double I1[N][N],I2[N][N];
    bzero(I1,sizeof(double)*N*N);
    bzero(I2,sizeof(double)*N*N);
        for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                I1[i][j]+=QT[i][k]*QT[j][k]; //QT*Q
                I2[i][j]+=QT[k][i]*QT[k][j]; //Q*QT
            }
        }
    }
    printf("I1=\n"); matprint(N,N,I1);
    printf("I2=\n"); matprint(N,N,I2);
    return 0;
}
