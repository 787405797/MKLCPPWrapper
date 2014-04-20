#include<iostream>
#include "util.h"
#include "math.h"
using namespace std;
#define ROW 2
#define COL 3
#define COL2 5

int main(){
    double x[ROW*COL]={1,1,1,0.5,0.5,0.5};
    double y[COL]={1,2,3};
    double z[ROW]={0,0};
    MKL_Complex16 cy[COL];
    MKL_Complex16 cx[ROW*COL];
    MKL_Complex16 cz[ROW];
    copy_VV(ROW,z,cz);
    copy_VV(COL,y,cy);
    for(int i=0;i<ROW*COL;i++){
        cx[i].real = cos(x[i]);
        cx[i].imag = sin(x[i]);
    }
    multiply_MV(ROW,COL,x,y,1,z);
    displayVec(ROW,1,z);

    displayVec(COL,ROW,cx);
    multiply_MtV(COL,ROW,cx,cy,1,cz);
    displayVec(ROW,1,cz);

    double f[8] = {1,2,3,4,5,6,7,8};
    MKL_Complex16 cf[8],cf_fft[16],cf_fft_ifft[16];
    for(int i=0;i<8;i++){
        cf[i].real = cos(f[i]);
        cf[i].imag = sin(f[i]);
    }
    cout<<"cf:"<<endl;
    displayVec(4,2,cf);
    if(fft(8,16,cf,cf_fft)){
        cout<<"Error!"<<endl;
        return 1;
    }
    cout<<"cf_fft:"<<endl;
    displayVec(4,4,cf_fft);

    if(ifft(16,16,cf_fft,cf_fft_ifft)){
        cout<<"Error!"<<endl;
        return 1;
    }
    cout<<"cf_fft_ifft:"<<endl;
    displayVec(4,4,cf_fft_ifft);

    cout<<"nextpow2(2999) = "<<nextpow2(2999)<<endl;
    cout<<"nextpow2(-0.1) = "<<nextpow2(-0.1)<<endl;
    cout<<"nextpow2(1e308) = "<<nextpow2(1e308)<<endl;
    cout<<"2<<11 = "<<(2<<11)<<endl;

    cout<<"cf:"<<endl;
    displayVec(1,8,cf);
    fftshift(8,cf);
    cout<<"cf (after shift):"<<endl;
    displayVec(1,8,cf);
    double f1[9] = {1,2,3,4,5,6,7,8,9};
    MKL_Complex16 cf1[9];
    for(int i=0;i<9;i++){
        cf1[i].real = cos(f[i]);
        cf1[i].imag = sin(f[i]);
    }
    cout<<"cf1:"<<endl;
    displayVec(1,9,cf1);
    fftshift(9,cf1);
    cout<<"cf1 (after shift):"<<endl;
    displayVec(1,9,cf1);
    
    return 0;
}
