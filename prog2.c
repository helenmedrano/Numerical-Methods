#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include <strings.h>

#define N_1 3
#define N_2 7

void matprint(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      printf("%g ",A[i][j]);
    }
    printf("\n");
  }
}

void LU(int size, double D[size][size],int n)
{
  int i,j,k,m,an,am;
  double x;

  printf("The matrix \n");
  matprint(size, size, D);

  for(k=0;k<n;k++)
  {
    for(j=k+1;j<n+1;j++)
    {
      x=D[j][k]/D[k][k];
      for(i=k;i<n+1;i++)
      {
        D[j][i]=D[j][i]-x*D[k][i];
      }
      D[j][k]=x;
    }
  }

  printf(" \n");
  printf("The matrix LU decomposed \n");
  matprint(size, size, D);
}

double vecnorm2(int n,double x[n]){
  double r=0;
  for(int i=0;i<n;i++){
    double t=fabs(x[i]);
    r+=t*t;
  }
  return sqrt(r);
}

void multAx(int m,int n, double A[m][n],double X[n],double B[m]){
  bzero(B,sizeof(double)*m);
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      B[i]+=A[i][j]*X[j];
    }
  }
}

double matnorm2(int n,double A[n][n]){
  double B[n][n],x[n],y[n],yk[n];
  bzero(B,sizeof(double)*n*n);
  for(int k=0;k<n;k++){
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        B[i][j]+=A[k][i]*A[k][j];
      }
    }
  }
  for(int k=0;k<n;k++){
    x[k]=2.0*random()/RAND_MAX-1.0;
  }
  multAx(n,n,B,x,y);
  double t,tk;
  for(int k=2;k<10;k++){
    multAx(n,n,B,y,yk);  // yk=B^k x
    t=vecnorm2(n,y);
    tk=vecnorm2(n,yk);
    //printf("%g,%g tk/t=%g\n",t,tk,tk/t);
    memcpy(y,yk,sizeof(double)*n);
  }
  return(sqrt(tk/t));
}

  /*  INVERSE */

  /* to find the inverse we solve [D][y]=[d] with only one element in 
  the [d] array put equal to one at a time */

double invmatnorm2(int size, int n, double D[size][size], double s[size][size]) {
  double d[size], y[size];
  double x;
  int m, i, j;

  for(m=0;m< size;m++)
  {
    d[0]=0.0;d[1]=0.0;d[2]=0.0;
    d[m]=1.0;
    for(i=0;i<=n;i++)
    {
      x=0.0;
      for(j=0;j<=i-1;j++)
      {
        x=x+D[i][j]*y[j];
      }
      y[i]=(d[i]-x);
    }

    for(i=n;i>=0;i--)
    {
      x=0.0;
      for(j=i+1;j<=n;j++)
      {
        x=x+D[i][j]*s[j][m];
      }
      s[i][m]=(y[i]-x)/D[i][i];
    }
  }

  /* Print the inverse matrix */
  printf("The Inverse Matrix\n");
  matprint(size, size, s);

  return matnorm2(size, s);
}

double matcond2(int size, double invmatnorm, double A[size][size]) {
  return invmatnorm*matnorm2(size, A);
}

int main(){
  int i,j,n_1,n_2,m,an,am;
  double a[N_2], d[N_1], B[N_2][N_2], C[N_1][N_1];
  double x,s_1[N_1][N_1],s_2[N_2][N_2],y_1[N_1], y_2[N_2];
  an=3; am=3;
  n_1 =2;
  n_2 = 6;

  printf("\nProblem solving for practice matrix:\n");
  /* Matrix 1*/
  double D[N_1][N_1] = {{7, 1, 7},
                        {5, 7, 7},
                        {3, 7, 8}};

  /* Copy of the original matrix 1*/
  memcpy(C, D, sizeof(double)*N_1*N_1);

  LU(N_1, D,n_1);

  double invmatnorm_practice = invmatnorm2(N_1, n_1, D, s_1);

  printf("Invmatnorm2 for practice problem: %g\n", invmatnorm_practice);
  printf("Cond for practice problem: %g\n", matcond2(N_1, invmatnorm_practice, C));

  printf("\nProblem solving for assigned matrix:\n");
  /* Matrix 2*/
  double A[N_2][N_2] = {
      { -6, -4,  6,  6, -2,  4,  7 },
      { -3, -4, -1,  8,  1, -2,  4 },
      { -1, -1, -5,  8, -6, -8, -7 },
      {  5,  2,  8, -9, -7, -4,  4 },
      {  4,  0, -5,  7,  6, -9,  3 },
      { -6,  5,  1,  0, -9, -9, -1 },
      {  2, -1, -6, -9,  7, -2,  8 }};

  double matnorm_special = matnorm2(N_2, A);


  /* Copy of the original matrix 1 in A to copy into B*/
  memcpy(B, A, sizeof(double)*N_2*N_2);

  LU(N_2, A,n_2);

  double invmatnorm_special = invmatnorm2(N_2, n_2, A, s_2);

  printf("matnorm2 for assigned matrix: %g\n", matnorm_special);
  printf("Invmatnorm2 for assigned matrix: %g\n", invmatnorm_special);
  printf("Cond for assigned matrix: %g\n", matcond2(N_2, invmatnorm_special, B));
}
