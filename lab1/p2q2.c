#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

double q6matnorm2(int n,double A[n][n]){
    double B[n][n],y[n],yk[n];
    bzero(B,sizeof(double)*n*n);
    for(int k=0;k<n;k++){                 // B = AT A
        for(int i=0;i<n;i++){   
            for(int j=0;j<n;j++){
                B[i][j]+=A[k][i]*A[k][j];
            }
        }
    }
    for(int i=0;i<n;i++){                 // Choose x in R^n randomly
        y[i]=2.0*random()/RAND_MAX+1.0;   // and store x in y for now
    }
    double q=0,qk;
    for(int k=1;k<100*n;k++){
        multAx(n,n,B,y,yk);               // yk = B^k x/|B^(k-1)x|
        qk=vecnorm2(n,yk);
        for(int j=0;j<n;j++){
            y[j]=yk[j]/qk;                // Overwrite y by y_k/|y_k|
        }
        if(fabs(qk-q)<5e-15*qk){          // Converge to 15 digits where
            return sqrt(qk);              // |A|=sqrt(|B^k x|/|B^(k-1)x|)
        }
        q=qk;
    }
    fprintf(stderr,"q6matnorm2: Failed to converge!\n");
    return sqrt(qk);
}

complex cq6matnorm2(int n,complex A[n][n]){
    complex B[n][n],y[n],yk[n];
    bzero(B,sizeof(complex)*n*n);
    for(int k=0;k<n;k++){                 // B = AT A
        for(int i=0;i<n;i++){   
            for(int j=0;j<n;j++){
                B[i][j]+=A[k][i]*A[k][j];
            }
        }
    }
    for(int i=0;i<n;i++){                 // Choose x in R^n randomly
        y[i]=2.0*random()/RAND_MAX+1.0;   // and store x in y for now
    }
    complex q=0,qk;
    for(int k=1;k<100*n;k++){
        cmultAx(n,n,B,y,yk);               // yk = B^k x/|B^(k-1)x|
        qk=cvecnorm2(n,yk);
        for(int j=0;j<n;j++){
            y[j]=yk[j]/qk;                // Overwrite y by y_k/|y_k|
        }
        
        if(cabs(qk-q)<cabs(5e-15*qk)){          // Converge to 15 digits where
            return sqrt(qk);              // |A|=sqrt(|B^k x|/|B^(k-1)x|)
        }
        q=qk;
    }
    fprintf(stderr,"q6matnorm2: Failed to converge!\n");
    return sqrt(qk);
}

double q6invmatnorm2(int n,double A[n][n]){
    double B[n][n],*P[n],y[n],yk[n];
    bzero(B,sizeof(double)*n*n);
    for(int k=0;k<n;k++){                 // B = AT A
        for(int i=0;i<n;i++){   
            for(int j=0;j<n;j++){
                B[i][j]+=A[k][i]*A[k][j];
            }
        }
    }
    PLUfact(n,B,P);
    for(int i=0;i<n;i++){                 // Choose x in R^n randomly
        y[i]=2.0*random()/RAND_MAX+1.0;   // and store x in y for now
    }
    double q=0,qk;
    for(int k=1;k<100*n;k++){
        PLUsolve(n,B,P,yk,y);             // yk = B^-k x/|B^(1-k)x|
        qk=vecnorm2(n,yk);
        for(int j=0;j<n;j++){
            y[j]=yk[j]/qk;                // Overwrite y by y_k/|y_k|
        }
        if(fabs(qk-q)<5e-15*qk){          // Converge to 15 digits where
            return sqrt(qk);              // |A|=sqrt(|B^-k x|/|B^(1-k)x|)
        }
        q=qk;
    }
    fprintf(stderr,"q6invmatnorm2: Failed to converge!\n");
    return sqrt(qk);
}

complex shiftinvpower(int n,complex A[n][n], complex alpha){
    double complex B[n][n],*P[n],y[n],yk[n];
    bzero(B,sizeof(complex)*n*n);
    
    /*for(int k=0;k<n;k++){                 // B = AT A
        for(int i=0;i<n;i++){   
            for(int j=0;j<n;j++){
                B[i][j]+=A[k][i]*A[k][j];
            }
        }
    }*/
    
    
    
    //PLUfact(n,B,P);
    cPLUfact(n, B, P);
    
    for(int i=0;i<n;i++){                 // Choose x in R^n randomly
        y[i]=2.0*random()/RAND_MAX+1.0;   // and store x in y for now
    }
    double complex q=0,qk;
    for(int k=1;k<100*n;k++){
        //PLUsolve(n,B,P,yk,y);             // yk = B^-k x/|B^(1-k)x|
        cPLUsolve(n, B, P, yk, y);
        
        qk=cvecnorm2(n,yk);
        for(int j=0;j<n;j++){
            y[j]=yk[j]/qk;                // Overwrite y by y_k/|y_k|
        }
        
        if(cabs(qk-q)<cabs(5e-15*qk)){          // Converge to 15 digits where
            return sqrt(qk);              // |A|=sqrt(|B^-k x|/|B^(1-k)x|)
        }
        q=qk;
    }
    fprintf(stderr,"q6invmatnorm2: Failed to converge!\n");
    return sqrt(qk);
}

double q6cond2(int n,double A[n][n]){
    return q6matnorm2(n,A)*q6invmatnorm2(n,A);
}

#define N 5
#include "matrix.c"
double X[N];

int main(){
    printf("Midterm Part 2 Question 2.\n");
    printf("N=%d\n",N);
    printf("A=\n"); cmatprint(N,N,A);
    for(int i=0;i<N;i++) X[i]=1;
    double t1,t2;  // for efficiency no need to compute norms twice.
    printf("|A|_2=%.15g\n",t1=cq6matnorm2(N,A));
    complex alpha = .5+4i;
    printf("|A^-1|_2=%.10g\n",t2=shiftinvpower(N,A, alpha));
//    printf("cond2(A)=%.10g\n",q6cond2(N,A));
    printf("cond2(A)=%.10g\n",t1*t2);
    return 0;
}
