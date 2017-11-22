#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

// Algorithm (a) from page 13 of Numerical Algorithms by Justin Solomon
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

  printf("Matrix LU:\n");
  matPrint(n, n, A);
  printf("Matrix P:\n");
  matPrint(n, n, P);
}

double vecnorm2(int n,double x[n]){
  double r=0;
  for(int i=0;i<n;i++){
    r+=x[i]*x[i];
  }
  return sqrt(r);
}

double matnorm2(int n,double A[n][n]){
  double x[n], XX[n];
  for(int i = 1; i <= n; i++) {
    x[i] = i;
  }

  double B[n][n],y[n],yk[n];
  bzero(B,sizeof(double)*n*n);
  for(int k=0;k<n;k++){ // B = ATA
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        B[i][j]+=A[k][i]*A[k][j];
      }
    }
  }
  for(int i=0;i<n;i++){ // Choose x ∈ Rn randomly
    y[i]=2.0*random()/RAND_MAX+1.0; // and store x in y for now
  }
  double q=0,qk;
  for(int k=1;k<100*n;k++){
    multAx(n,n,B,y,yk); // yk = Bk
    qk=vecnorm2(n,yk);
    for(int j=0;j<n;j++){
      y[j]=yk[j]/qk; // Overwrite y by yk/∥yk∥2
    }
    if(fabs(qk-q)<5e-15*qk){ // Converge to 15 digits where
      return sqrt(qk); // ∥A∥2 ≈ (∥Bk
    }
    q=qk;
  }
  fprintf(stderr,"matnorm2: Failed to converge!\n");
  return sqrt(qk);
}

void PLUinverse(int n, double LU[n][n], double P[n][n], double s[n][n]) {
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      if (P[i][j] == 1)
        s[i][j] = 1.0;
      else
        s[i][j] = 0.0;

      for (int k = 0; k < i; k++)
        s[i][j] -= LU[i][k] * s[k][j];
    }

    for (int i = n - 1; i >= 0; i--) {
      for (int k = i + 1; k < n; k++) {
        s[i][j] -= LU[i][k] * s[k][j];
      }
      s[i][j] = s[i][j] / LU[i][i];
    }
  }
  printf("\ninverse: \n");
  matPrint(n, n, s);

}

#define M 3
#define N 2
#define Q 4

int main() {
  double A[M][N] = {{1, 2},
                    {3, 4},
                    {5, 6}};
  double B[N][Q] = {{1, 2, 3, 4},
                    {5, 6, 7, 8}};
  double C[M][M] = {{7, 1, 7},
                    {5, 7, 7},
                    {3, 7, 8}};
  double A_mult_B[M][Q];
  double s[M][M];

  printf("Matrix A:\n");
  matPrint(M, N, A);

  printf("Matrix B:\n");
  matPrint(N, Q, B);

  matrixMult(M, N, Q, A, B, A_mult_B);
  printf("A_mult_B:\n");
  matPrint(M, Q, A_mult_B);

  printf("Matrix C:\n");
  matPrint(M, M, C);

  /* PLU factorization on 3x3 matrix */
  double P[M][M];
  PLUfact(M, C, P);

  PLUinverse(M, C, P, s);
  printf("matnorm2 for inverse of matrix: %g\n", matnorm2(M, s));

  return 0;
}
