#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <strings.h>
#include <cilk/cilk.h>
#include <lapacke.h>
#include "matrixlib.h"

double dotprod(int n,double x[n],double y[n]){
    double s=0;
    for(int i=0;i<n;i++){
        s+=x[i]*y[i];
    }
    return s;
}

double phi0(double x){ return 1.0; }
double phi1(double x){ return x; }
double phi2(double x){ return x*x; }
double phi3(double x){ return x*x*x; }

typedef double func(double);
func *phi[4] = { phi0, phi1, phi2, phi3 };

void gramschmidt(int m,int n,double A[m][n],
    double QT[n][m]){
    double vtilde[m];
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++){
            vtilde[i]=A[i][j];
        }
        for(int k=0;k<j;k++){
            double adotv=dotprod(m,vtilde,QT[k]);
            for(int l=0;l<m;l++){
                vtilde[l]-=adotv*QT[k][l];
            }
        }
        double vnorm=vecnorm2(m,vtilde);
        for(int i=0;i<m;i++){
            QT[j][i]=vtilde[i]/vnorm;
        }
    }
}

void hheliminate(int m,int n,double w[m],
    double A[m][n]){
    for(int j=0;j<n;j++){
        double at[m];
        for(int i=0;i<m;i++) at[i]=A[i][j];
        double wdpat=dotprod(m,w,at);
        for(int i=0;i<m;i++) at[i]-=2*w[i]*wdpat;
        for(int i=0;i<m;i++) A[i][j]=at[i];
    }
}
void householder(int m,int n,double A[m][n],
    double R[m][n], double Q[m][n]){
    double Hs[n][m];
    memcpy(R,A,sizeof(double)*m*n);
    for(int j=0;j<n;j++){
        double vtilde[m];
        for(int i=0;i<m;i++){
            if(i<j) vtilde[i]=0;
            else vtilde[i]=R[i][j];
        }
        double c=vecnorm2(m,vtilde);
        if(vtilde[j]<0) vtilde[j]-=c;
        else vtilde[j]+=c;
        double vnorm=vecnorm2(m,vtilde);
        for(int i=0;i<m;i++){
            vtilde[i]/=vnorm;
        }
        hheliminate(m,n,vtilde,R);
 //       printf("R=\n"); matprint(m,n,R);
        memcpy(Hs[j],vtilde,sizeof(double)*m);
    }
    bzero(Q,sizeof(double)*m*n);
    for(int j=0;j<n;j++) Q[j][j]=1;
    for(int j=n-1;j>=0;j--){
        hheliminate(m,n,Hs[j],Q);
    }
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            if(i>j) R[i][j]=0;
        }
    }
}

#define N 4

int main(){
    setrlimit(RLIMIT_STACK,
        &(const struct rlimit)
        {RLIM_INFINITY,RLIM_INFINITY});
    printf("N=%d\n",N);
    FILE *fp=fopen("file05a.dat","r");
    int m;
    fscanf(fp,"%d",&m);
    printf("m=%d\n",m);
    double A[m][N],Q[m][N],R[m][N];
    double X[m],Y[m];
    for(int i=0;i<m;i++){
        fscanf(fp,"%lf %lf",&X[i],&Y[i]);
        printf("%lf %lf\n",X[i],Y[i]);
    }
    for(int i=0;i<m;i++){
        for(int j=0;j<N;j++){
            A[i][j]=phi[j](X[i]);
        }
    }
    householder(m,N,A,R,Q);
    double qty[N];
    bzero(qty,sizeof(double)*N);
    for(int i=0;i<m;i++){
        for(int j=0;j<N;j++){
            qty[j]+=Q[i][j]*Y[i];
        }
    }
    double c[N];
    backsub(N,R,c,qty);
    printf("c=\n"); vecprint(N,c);
    return 0;
}    
