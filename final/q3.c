#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>

#define N 4

/*
 * w is above the diagonal == c
 * v is on the diagonal == b
 * u is below diagonal == a
 * solve for x
 * b is given == d
*/
void tdthomas(int n,double u[n],double v[n],double w[n-1], double x[n],double b[n]) {
  memcpy(x, b, sizeof(double)*n);
  n--;
  w[0] /= v[0];
  x[0] /= v[0];

  for (int i = 1; i < n; i++) {
    w[i] /= v[i] - u[i]*w[i-1];
    x[i] = (x[i] - u[i]*x[i-1]) / (v[i] - u[i]*w[i-1]);
  }

  x[n] = (x[n] - u[n]*x[n-1]) / (v[n] - u[n]*w[n-1]);

  for (int i = n; i-- > 0;) {
    x[i] -= w[i]*x[i+1];
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

int main() {
  double u[N] = { 0, -1, -1, -1 };
  double v[N] = { 4,  4,  4,  4 };
  double w[N-1] = {-1, -1, -1};
  double b[N] = { 5,  5, 10, 23 };
  double x[N];
  // results    { 2,  3,  5, 7  }
  tdthomas(N, u, v, w, x, b);
  vecprint(N, x);

  return 0;
}
