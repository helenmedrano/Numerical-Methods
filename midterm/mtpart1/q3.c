/*  Programming template for Midterm Part 1 Question 3

    Please complete the four routines matmatAB, matmatATB, matmatABT
    and matmatATBT to multiply the matrices A and B and their transposes
    such that the four resulting tests in main produce the messages
    "Correct!  matmatXXX passed the rounding error test!" in every case.

    You may add additional code for debugging and to modularize your
    computation; however, there is no need to do so.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <strings.h>
#include <math.h>

void transposeMatrix(int m, int n, double A[m][n], double result[n][m]) {
  for(int i = 0; i < m; i++) {
      for(int j = 0; j < n; j++) {
            result[j][i] = A[i][j];
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

void matPrint(int m, int n, double A[m][n]) {
  for(int i = 0; i < m; i++) {
      for(int j = 0; j < n; j++) {
            printf("%g ", A[i][j]);
          }
      printf("\n");
    }
}

//  Compute the distance ||A-B||_frob between to matrices.
double distfrob(int m,int n,double A[m][n],double B[m][n]){
    double r=0;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        double t=A[i][j]-B[i][j];
        r+=t*t;
    }
    return sqrt(r);
}

//  Multiply the matrix A by B and return the result as AB.
void matmatAB(int m,int p,int n,double A[m][p],double B[p][n],
    double AB[m][n]){
    matrixMult(m, p, n, A, B, AB);

/*  Fill in the code here. */

}

//  Multiply the matrix AT by B and return the result as ATB.
void matmatATB(int m,int p,int n,double A[p][m],double B[p][n],
    double ATB[m][n]){
    double A_T[m][p];
    transposeMatrix(p, m, A, A_T);
    matrixMult(m, p, n, A_T, B, ATB);

/*  Fill in the code here. */

}

//  Multiply the matrix A by BT and return the result as ABT.
void matmatABT(int m,int p,int n,double A[m][p],double B[n][p],
    double ABT[m][n]){
    double B_T[p][n];
    transposeMatrix(n, p, B, B_T);
    matrixMult(m, p, n, A, B_T, ABT);

/*  Fill in the code here. */

}

//  Multiply the matrix AT by BT and return the result as ATBT.
void matmatATBT(int m,int p,int n,double A[p][m],double B[n][p],
    double ATBT[m][n]){
    double A_T[m][p];
    transposeMatrix(p, m, A, A_T);

    double B_T[p][n];
    transposeMatrix(n, p, B, B_T);

    matrixMult(m, p, n, A_T, B_T, ATBT);

/*  Fill in the code here. */

}

double A[4][5]={
    { -2,  2,  1, -1,  2 },
    {  2,  2, -2,  2, -2 },
    { -1, -2,  2,  1, -2 },
    {  1, -1, -0,  0, -0 }};

double B[5][6]={
    { -2, -0, -1, -2,  2, -2 },
    {  0, -2, -2, -2,  0,  1 },
    {  2, -1, -1,  2, -2, -2 },
    {  1, -2, -0, -2,  0, -1 },
    { -1, -2, -1, -2, -2,  1 }};

double C[4][6]={
    {  2,  0, -1,  2, -0,  2 },
    {  0,  1, -2, -2,  2,  2 },
    {  0, -2,  2, -0,  1, -1 },
    {  2, -1, -2,  1, -2, -1 }};

double AB[4][6]={
    {  3,  -7,  -5,   0, -10,   7 },
    { -4,  -2,  -2, -12,  12,  -2 },
    {  9,   4,   5,  12,  -2,  -7 },
    { -2,   2,   1,   0,   2,  -3 }};

double CTA[6][5]={
    { -2,   2,   2,  -2,   4 },
    {  3,   7,  -6,   0,   2 },
    { -6,  -8,   7,  -1,  -2 },
    { -7,  -1,   6,  -6,   8 },
    {  1,   4,  -2,   5,  -6 },
    {  0,  11,  -4,   1,   2 }};

double CBT[4][5]={
    { -11,   0,   5,  -4,  -3 },
    {   6,   8, -11,   0,   2 },
    {   2,  -1,   0,   5,  -1 },
    {  -6,   3,  15,   3,   3 }};

double BTAT[6][4]={
    {   3,  -4,   9,  -2 },
    {  -7,  -2,   4,   2 },
    {  -5,  -2,   5,   1 },
    {   0, -12,  12,   0 },
    { -10,  12,  -2,   2 },
    {   7,  -2,  -7,  -3 }};

double zAB[4][6],zCTA[6][5],zCBT[4][5],zBTAT[6][4];

int main(){
    printf("Midterm Part 1 Question 3.\n");
    double rho;
    matmatAB(4,5,6,A,B,zAB);
    rho=distfrob(4,6,AB,zAB);
    printf("||AB-matmat|| = %.14e\n",rho);
    if(rho<1e-7){
        printf("Correct!  matmatAB passed the rounding error test!\n");
    } else {
        printf("Incorrect!  matmabAB did not pass the rounding test!\n");
    }
    matmatATB(6,4,5,C,A,zCTA);
    rho=distfrob(6,5,CTA,zCTA);
    printf("||CTA-matmat|| = %.14e\n",rho);
    if(rho<1e-7){
        printf("Correct!  matmatATB passed the rounding error test!\n");
    } else {
        printf("Incorrect!  matmabATB did not pass the rounding test!\n");
    }
    matmatABT(4,6,5,C,B,zCBT);
    rho=distfrob(4,5,CBT,zCBT);
    printf("||CBT-matmat|| = %.14e\n",rho);
    if(rho<1e-7){
        printf("Correct!  matmatABT passed the rounding error test!\n");
    } else {
        printf("Incorrect!  matmabABT did not pass the rounding test!\n");
    }
    matmatATBT(6,5,4,B,A,zBTAT);
    rho=distfrob(6,4,BTAT,zBTAT);
    printf("||BTAT-matmat|| = %.14e\n",rho);
    if(rho<1e-7){
        printf("Correct!  matmatATBT passed the rounding error test!\n");
    } else {
        printf("Incorrect!  matmabATBT did not pass the rounding test!\n");
    }
    return 0;
}
