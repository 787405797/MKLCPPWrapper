#include<iostream>
#include<iomanip>
#include <cmath>
#include "util.h"
#include <mkl.h>
using namespace std;

const CBLAS_ORDER order = CblasRowMajor;
const CBLAS_TRANSPOSE NoTrans = CblasNoTrans;
const CBLAS_TRANSPOSE Trans = CblasTrans;
const float beta = 0;
const MKL_INT incX = 1;
const MKL_INT incY = 1;

// VecOut(Row) = Mat(Row*Col) * Vec(Col) * alpha
void multiply_MV(int Row, int Col, const double *Mat, const double *Vec,const float alpha, double *VecOut){
    cblas_dgemv(order, NoTrans, Row, Col, alpha, Mat, Col, Vec, incX, beta, VecOut, incY );
}

void multiply_MtV(int Row, int Col,const MKL_Complex16 *Mat,const MKL_Complex16 *Vec, const float alpha,MKL_Complex16 *VecOut){
    MKL_Complex16 calpha,cbeta;
    calpha.real = alpha;
    calpha.imag = 0;
    cbeta.real = beta;
    cbeta.imag = 0;
    cblas_zgemv(order,Trans, Row, Col, &calpha, Mat, Col, Vec, incX, &cbeta, VecOut, incY);
}
// VecOut = vec1 + vec2
// vec2 = vec1 + vec2
void add_VV(int N,const double *vec1, double * vec2, double *VecOut=NULL){
    if (VecOut == NULL){
        double *VecOut = new double[N];
        vdAdd(N, vec1, vec2, VecOut);
        copy_VV(N, VecOut, vec2);
        delete [] VecOut;
    }
    else{
        vdAdd(N, vec1, vec2, VecOut);
    }
}

void copy_VV(int N, const double *VecIn, double *VecOut){
    cblas_dcopy(N, VecIn, incX, VecOut, incY);
}

void copy_VV(int N,const double *VecIn,MKL_Complex16 *VecOut){
    for(int i=0;i<N;i++){
        VecOut[i].real = VecIn[i];
        VecOut[i].imag = 0;
    }
}

void copy_VV(int N,const MKL_Complex16 *VecIn,MKL_Complex16 *VecOut){
    for(int i=0;i<N;i++){
        VecOut[i].real = VecIn[i].real;
        VecOut[i].imag = VecIn[i].imag;
    }
}

void displayVec(int Row, int Col, const double *Vec){
    for(int i=0;i<Row;i++){
        for(int j=0;j<Col;j++){
            cout<<Vec[i*Col + j]<<"\t";            
        }
        cout<<endl;
    }
}

void displayVec(int Row, int Col, const MKL_Complex16 *Vec){
    for(int i=0;i<Row;i++){
        for(int j=0;j<Col;j++){
            cout<<resetiosflags(ios::showpos)<<Vec[i*Col + j].real<<setiosflags(ios::showpos)<<Vec[i*Col+j].imag<<"j"<<"\t";
        }
        cout<<endl;
    }
    
}

int fft(int N, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut){
    copy_VV(N,VecIn,VecOut);
    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX,1,(MKL_LONG)N);
    if (status != 0) return status;
    status = DftiCommitDescriptor(hand);
    if (status != 0) return status;
    status = DftiComputeForward(hand,VecOut);
    if (status != 0) return status;
    DftiFreeDescriptor(&hand);
    return status;
}

int fft(int N, int fftNum, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut){
    if (fftNum <= N)
        return fft(fftNum,VecIn,VecOut);
    MKL_Complex16 *temp = new MKL_Complex16[fftNum];
    copy_VV(N,VecIn,temp);
    for(int i=N;i<fftNum;i++){
        temp[i].real = 0;
        temp[i].imag = 0;
    }
    return fft(fftNum,temp,VecOut);
}

int ifft(int N, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut){
    copy_VV(N,VecIn,VecOut);
    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX,1,(MKL_LONG)N);
    if (status != 0) return status;
    status = DftiCommitDescriptor(hand);
    if (status != 0) return status;
    status = DftiComputeBackward(hand,VecOut);
    cblas_zdscal(N,1/(double)N,VecOut,1);
    if (status != 0) return status;
    DftiFreeDescriptor(&hand);
    return status;
}

int ifft(int N, int fftNum, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut){
    if (fftNum <= N)
        return ifft(fftNum,VecIn,VecOut);
    MKL_Complex16 *temp = new MKL_Complex16[fftNum];
    copy_VV(N,VecIn,temp);
    for(int i=N;i<fftNum;i++){
        temp[i].real = 0;
        temp[i].imag = 0;
    }
    return ifft(fftNum,temp,VecOut);
}

void fftshift(int N,MKL_Complex16 *VecInOut){
    if(N%2){
        MKL_Complex16 *temp = new MKL_Complex16[N];
        copy_VV(N/2+1,VecInOut,temp+N/2);
        copy_VV(N/2,VecInOut+N/2+1,temp);
        copy_VV(N,temp,VecInOut);
    }
    else{
        cblas_zswap(N/2,VecInOut,1,VecInOut+N/2,1);
    }
}

int nextpow2(const double x){
    // Assuming that x is bigger than 0;
    return ceil(log10(x)/log10(2));
}

