#ifndef __UTIL_H
#define __UTIL_H
#include "mkl.h"

// VecOut(Row) = Mat(Row*Col) * Vec(Col)
void multiply_MV(int Row, int Col, const double *Mat, const double *Vec, const float alpha, double *VecOut);
void multiply_MtV(int Row, int Col, const MKL_Complex16 *Mat, const MKL_Complex16 *Vec, const float alpha, MKL_Complex16 *VecOut);


// VecOut = vec1 + vec2
// vec2 = vec1 + vec2
void add_VV(int N,double *vec1, double * vec2, double *VecOut=NULL);

// VecOut = VecIn
void copy_VV(int N, const double *VecIn, double *VecOut);
void copy_VV(int N, const double *VecIn,MKL_Complex16 *VecOut);
void copy_VV(int N, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut);

// Display vector or matrix
void displayVec(int Row, int Col, const double *Vec);
void displayVec(int Row, int Col, const MKL_Complex16 *Vec);

// fft and ifft
int fft(int N, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut);
int fft(int N, int fftNum, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut);
int ifft(int N, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut);
int ifft(int N, int fftNum, const MKL_Complex16 *VecIn, MKL_Complex16 *VecOut);
void fftshift(int N,MKL_Complex16 *VecInOut);

// other Matlab function
int nextpow2(const double x);

#endif
