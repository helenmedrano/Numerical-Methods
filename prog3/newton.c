#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 9

// x^3-3
double function(double x) {
  return x*x*x-3;
}

// 3x^2
double derivative(double x) {
  return 3*x*x;
}

// e_n = x_n - 3^1/3
double e_n(double x) {
  return x-pow(3, 1/3);
}

void printTableRow(int n, double x_n, double e_n) {
  printf("%d | %lf | %lf\n", n, x_n, e_n);
}

void newtonsMethod(double x, double e[N]) {
  double h, x_n;

  printf("n | x_n | e_n\n");

  for(int i = 0; i < N; i++) {
    h = function(x)/derivative(x);
    x_n = x-h;

    e[i] = e_n(x_n);
    printTableRow(i, x_n, e[i]);

    x = x_n;
  }

  printf("The root's value is: %f\n", x);
}

// |e_n+1| = M_n*|e_n|^2
// M_n = |e_n+1|/|e_n|^2
void M_n(double e_n[N], double M[N]) {
  int n = 4;
  printf("n | M_n\n");
  for(int i = 1; i <= n; i++) {
    M[i] = fabs(e_n[i+1])/(pow(fabs(e_n[i]), 2));
    printf("%d | %f\n", i, M[i]);
  }
}

int main() {
  double x_0 = 1;
  double n[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  double e_n[N], m_n[N];
  newtonsMethod(x_0, e_n);
  M_n(e_n, m_n);

  return 0;
}
