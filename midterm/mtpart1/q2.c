/*  Programming template for Midterm Part 1 Question 2

    Please complete the routine matinv below with code to compute
    the inverse of matrix A so that the resulting test in main 
    produces the message "Correct!  Passed the rounding error test!"

    You may use plufact and plusolve in your program but note that
    plufact overwrites the input A with the LU factorization.

    You may add additional code for debugging and to modularize your
    computation; however, there is no need to do so.

    If you change the matrix size from N=256 to something smaller
    for debugging, please change it back before submitting.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <strings.h>
#include <math.h>

//  Compute the Frobenius matrix norm.
double matfrobnorm(int n,double A[n][n]){
    double r=0;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++){
        double t=A[i][j];
        r+=t*t;
    }
    return sqrt(r);
}

//  Compute the PLU factorization of A overwriting A with L and U.
//  Upon exit U is stored in the diagonal of A and the upper triagular
//  part while L is stored as the lower triangular part of A with the
//  diagonal assume to be 1.  The permutation matrix P is given by a
//  vector of pointers.
void plufact(int n,double A[n][n],double *P[n]){
    for(int i=0;i<n;i++){
        P[i]=&A[i][0];
    } for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
            if(fabs(P[i][i])<fabs(P[j][i])){
                double *t=P[i]; P[i]=P[j]; P[j]=t;
            }
        } for(int j=i+1;j<n;j++){
            double alpha=P[j][i]/P[i][i]; for(int k=i+1;k<n;k++){
                P[j][k]-=alpha*P[i][k];
            } P[j][i]=alpha;
        }
    }
}

//  Use the PLU factorization found by plufact above to solve the linear
//  equation PLU x = b for x.
void plusolve(int n,double LU[n][n],double *P[n],
    double x[n],double b[n]){
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

//  Find the inverse of matrix A and return it as B.  This routine does
//  not overwrite or change the value of A.
void matinv(int n,double A[n][n],double B[n][n]){

/*  Please fill this routine with code to find the inverse of matrix A.
    You may use the plufact and plusolve routines above, but be careful
    not to overwrite the values of matrix A.  For example, you may want
    to copy A to a temporary variable using code such as

        double LU[n][n];
        memcpy(LU,A,sizeof(double)*n*n);

    before calling plufact.  */

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
    printf("\nMatrix A: \n");
    matPrint(n, n, A);

    printf("\ninverse: \n");
    matPrint(n, n, B);
}




//#define N 256
#define N 4
double A[N][N],Ainv[N][N],R[N][N];
int main(){
    printf("Midterm Part 1 Question 2.\n");
    for(int i=0;i<N;i++) for(int j=0;j<N;j++){
        A[i][j]=2.0*random()/RAND_MAX-1.0;
    }
    matinv(N,A,Ainv);
    //  At this point A is still the same as it was before.
    bzero(R,sizeof(double)*N*N);
    for(int i=0;i<N;i++){
        for(int k=0;k<N;k++){
            for(int j=0;j<N;j++){
                R[i][j]+=A[i][k]*Ainv[k][j];
            }
        }
    }
    //  At this point R should approximate the identity matrix
    for(int i=0;i<N;i++) R[i][i]-=1.0;
    //  At this point R should approximate zero matrix
    double rho=matfrobnorm(N,R);
    printf("|A*Ainv-I|_frob = %.14e\n",rho);
    if(rho<1e-7){
        printf("Correct!  Passed the rounding error test!\n");
    } else {
        printf("Incorrect!  Did not pass the rounding error test!\n");
    }
    return 0;
}
