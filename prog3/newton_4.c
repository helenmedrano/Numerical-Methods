#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 4

// x^3-3
double function(double x) {
  return x*x*x-3;
}

// 3x^2
double derivative(double x) {
  return 3*x*x;
}

double g_x(double x) {
  return x-function(x)/derivative(x);
}

double e_n(double x, double x_inf) {
  return x-x_inf;
}

double M_n(double e_nPlusOne, double e) {
  return fabs(e_nPlusOne/e/e);
}

int main() {
  double x_0 = 1;
  double x = x_0;
  double y, e, e_nPlusOne, M;
  double x_inf = pow(3., 1./3.);
  printf("%2s %22s %22s %22s\n","n","x_n","e_n", "M_n");
  for(int i = 0; i <= N ; i++) {
    e = e_n(x, x_inf);
    y = g_x(x);
    e_nPlusOne = e_n(y, x_inf);
    M = M_n(e_nPlusOne, e);
    printf("%3d %22.14e %22.14e %22.14e\n", i, x, e, M);
    x = y;
  }

  return 0;
}
