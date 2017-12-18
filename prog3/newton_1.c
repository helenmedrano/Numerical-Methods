#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 100000

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

int main() {
  double x_0 = 1;
  double x = x_0;
  double x_inf = pow(3, 1./3.);
  printf("xinf: %22.14e\n", x_inf);
  printf("%2s %22s %22s\n","n","xn","en");
  for(int i = 0; i < N ; i++) {
    printf("%3d %22.14e %22.14e\n", i, x, e_n(x, x_inf));
    x = g_x(x);
  }

  return 0;
}
