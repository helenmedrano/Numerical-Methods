/*  Programming template for Midterm Part 1 Question 1

    Please complete the routine ifft below with an inverse fast
    Fourier routine so that the resulting test in main produces the
    message "Correct!  Passed the rounding error test!"

    You may add additional code for debugging and to modularize your
    computation; however, there is no need to do so.

    If you change the transform size from N=256 to something smaller
    for debugging, please change it back before submitting.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

// Compute the distance between two complex vectors of length n.
double cdist2(int n,complex x[n],complex y[n]){
    double r=0;
    for(int k=0;k<n;k++) {
        double t=x[k]-y[k];
        r+=t*conj(t);
    }
    return sqrt(r);  // ||x-y||_2
}

// Compute the fast Fourier transform of a vector x of length n
// and return the result in b.  Call this routine with s=1.
void fft(int n,int s,complex x[],complex b[]){
    if(n==1){
        b[0]=x[0];
        return;
    }
    if(n%2){
        printf("Error: n was not a power of 2!\n");
        exit(1);
    }
    int K=n/2;
    fft(K,2*s,&x[s],&b[K]);
    fft(K,2*s,&x[0],&b[0]);
    for(int k=0;k<K;k++){
        complex w=cexp(-I*2*M_PI*k/n);
        complex be=b[k],bo=b[k+K];
        b[k]=be+w*bo;
        b[k+K]=be-w*bo;
    }
}

// Compute the inverse Fourier transform of a vector x of length n
// and return the result in b.  Note call this routine with s=1.
void ifft(int n,int s,complex x[],complex b[]){

/*  Please fill this routine with code to find the inverse discrete
    Fourier transform of a vector with length equal to a power of two
    using a fast conquer-and-divide algorithm. */
    if(n==1){
        b[0]=x[0];
        return;
    }
    if(n%2){
        printf("Error: n was not a power of 2!\n");
        exit(1);
    }
    int K=n/2;
    ifft(K,2*s,&x[s],&b[K]);
    ifft(K,2*s,&x[0],&b[0]);
    for(int k=0;k<K;k++){
        complex w=cexp(I*2*M_PI*k/n)/n;
        complex be=b[k],bo=b[k+K];
        b[k]=be+w*bo;
        b[k+K]=(be-w*bo);
    }


}

void cmatprint(int m,int n,complex A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      printf("(%g %g) ",A[i][j]);
    }
    printf("\n");
  }
}
void cvecprint(int m,complex X[m]){
  cmatprint(m,1,(complex (*)[1])X);
}

#define N 256
//#define N 4
complex X[N],Y[N],Z[N];
int main(){
    printf("Midterm Part 1 Question 1.\n");
    for(int k=0;k<N;k++) {
        X[k]=cexp(2*M_PI*I*random()/RAND_MAX)*random()/RAND_MAX;
    }
    fft(N,1,X,Y);
    ifft(N,1,Y,Z);
    //  At this point X and Z should be identical upto rounding error
    double rho=cdist2(N,X,Z);
    printf("|X-ifft(fft(X))|_2 = %.14e\n",rho);
    if(rho<1e-7){
        printf("Correct!  Passed the rounding error test!\n");
    } else {
        printf("Incorrect!  Did not pass the rounding error test!\n");
    }
    return 0;
}
