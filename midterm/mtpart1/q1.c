#include <stdio.h>
#include <strings.h>
#include <complex.h>
#include <math.h>
#define N 4
double A[N];
double X[N];

void matprint(int m,int n,double A[m][n]){
    for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                        printf("%g ",A[i][j]);
                    }
            printf("\n");
        }
}
void vecprint(int m,double X[m]){
    matprint(m,1,(double (*)[1])X);
}

void multiplyAx(int n,double A[n], double X[n],double B[n]){
    bzero(B,n*sizeof(double));
    for(int i=0;i<n;i++){
        B[i]+=A[i]*X[i];
    }
}

void initX() {
	for(int k=0;k<N;k++){
		X[k]=1.0/(1+k)*(2+k)+I*1.0/(N-k);
	}
}
void dft(int n, double A[n]) {
	for(int k=0;k<N;k++){
		A[k]=cexp(-I*2*M_PI*k/N);
	}
}
void pdft(int n, int s, double x[n], double A[n]){
	if(N%2) {
		printf("Error: n is not a power of 2!\n");
		exit(1);
	}
	if(n==1) {
		A[0] = x[0];
		return;
	}
	int K=n/2;
	pdft(K,2*s,&x[s],&A[K]);
	pdft(K,2*s,&x[0],&A[0]);
	for(int k=0;k<K;k++){
		double w=cexp(-I*2*M_PI*k/n);
		A[k]=w;
		A[k+K]=w;
	}
}
int main(){
	initX();
	double B[N];
	// DFT
	printf("============= DFT ==============\n");
    dft(N, A);
    // double C[N];
    printf("A=\n"); vecprint(N,A);
    //printf("AC=\n"); matprint(N,N,AC);
    printf("X=\n"); vecprint(N,X);
    multiplyAx(N,A,X,B);  // B=A*X
    printf("B=\n"); vecprint(N,B);
    
	// pdft
	printf("============= PDFT ==============\n");
	pdft(N, 1, X, A);
    // double C[N];
    printf("A=\n"); vecprint(N,A);
    //printf("AC=\n"); matprint(N,N,AC);
    printf("X=\n"); vecprint(N,X);
    multiplyAx(N,A,X,B);  // B=A*X
    printf("B=\n"); vecprint(N,B);
    return 0;
}
