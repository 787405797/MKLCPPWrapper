#include "util.h"
#include "LFM_radar.h"
#include <math.h>
#include <malloc.h>
#include <mkl.h>
#include <iostream>
#include <iomanip>
using namespace std;

void LFM_radar(double T, double B, double Rmin, double Rmax, double* R, double* RCS, int N){
    // parameters
    double C = 3e8; // propagation speed
    double K = B/T; // chirp slope
    double Rwid = Rmax - Rmin; // receive window in meter
    double Twid = 2 * Rwid / C; // receive window in second
    double Fs = 5 * B, Ts = 1 / Fs; // sampling frequency and sampling spacing
    int Nwid = ceil(Twid / Ts); // receive window in number

    // Generate the echo
    double *t = (double *)malloc(sizeof(double) * Nwid); // receive window
    for(int i=0; i<Nwid; i++){
        t[i] = 2*Rmin/C + i*(2*Rmax/C - 2*Rmin/C)/(Nwid - 1);
    }
    
    int M = N; // number of targets

    double *td = (double *)malloc(sizeof(double) * M * Nwid);
    for(int i=0;i<M;i++){
        for(int j=0;j<Nwid;j++){
            td[i*Nwid+j] = t[j] - 2*R[i]/C;
        }
    }

    MKL_Complex16 *temp = new MKL_Complex16[M*Nwid];
    for(int i=0;i<M*Nwid;i++){
        if(fabs(td[i]) < T/2){
            temp[i].real = cos(M_PI * K * td[i] * td[i]);
            temp[i].imag = sin(M_PI * K * td[i] * td[i]);
        }
        else{
            temp[i].real = 0;
            temp[i].imag = 0;
        }
    }
    
    MKL_Complex16 *cRCS = new MKL_Complex16[M];
    MKL_Complex16 *cSrt = new MKL_Complex16[Nwid];
    copy_VV(M,RCS,cRCS);
    multiply_MtV(M,Nwid,temp,cRCS,1,cSrt);    

    /*
     * Digital processing of pulse compression radar using FFT and IFFT
     */
    int Nchirp = ceil(T/Ts); // pulse duration in number
    long long Nfft = 2<<(nextpow2(Nwid+Nwid-1)-1);
    MKL_Complex16 * cSrw = new MKL_Complex16[Nfft]; // fft of radar echo
    fft(Nwid,Nfft,cSrt,cSrw); 
    double *t0 = new double[Nchirp];
    for(int i=0;i<Nchirp;i++)
        t0[i] = -T/2 + (T/(Nchirp-1))*i;
    MKL_Complex16 * cSt = new MKL_Complex16[Nchirp]; // Chirp signal
    for(int i=0;i<Nchirp;i++){
        cSt[i].real = cos(M_PI*K*t0[i]*t0[i]);
        cSt[i].imag = sin(M_PI*K*t0[i]*t0[i]);
    }
    MKL_Complex16 *cSw = new MKL_Complex16[Nfft];
    fft(Nchirp,Nfft,cSt,cSw);
    MKL_Complex16 *cSot = new MKL_Complex16[Nfft];
    MKL_Complex16 *ctemp1 = new MKL_Complex16[Nfft];
    vzMulByConj(Nfft,cSrw,cSw,ctemp1);
    
    cout<<"ctemp1:"<<endl;
    displayVec(20,5,ctemp1);
    cout<<"..."<<endl;

    ifft(Nfft,ctemp1,cSot);

    cout<<"cSot(before shift):"<<endl;
    displayVec(20,5,cSot);
    cout<<"..."<<endl;

    fftshift(Nfft,cSot);
    
    /*
     * print cSot and cSrt
     */
    cout<<setprecision(4);
    cout<<"td:"<<endl;
    displayVec(20,5,td);
    cout<<"..."<<endl;
    cout<<"cSrt:"<<endl;
    displayVec(20,5,cSrt);
    cout<<"..."<<endl;
    cout<<"Nfft:"<<endl;
    cout<<Nfft<<endl;
    cout<<"cSrw:"<<endl;
    displayVec(20,5,cSrw);
    cout<<"..."<<endl;
    cout<<"t0:"<<endl;
    displayVec(20,5,t0);
    cout<<"..."<<endl;
    cout<<"cSt:"<<endl;
    displayVec(20,5,cSt);
    cout<<"..."<<endl;
    cout<<"cSw:"<<endl;
    displayVec(20,5,cSw);
    cout<<"..."<<endl;
    cout<<"cSot:"<<endl;
    displayVec(20,5,cSot);
    cout<<"..."<<endl;
    
    /*
     * free source
     */
    free(cSot);
    free(ctemp1);
    free(cSw);
    free(cSt);
    free(cSrw);
    free(cRCS);
    free(cSrt);
    free(temp);
    free(t);
}

void LFM_radar(){
    double T = 10e-6; // 10us
    double B = 30e6;  // chirp frequency modulation bandwidth 30MHz
    double Rmin = 10000, Rmax = 15000; // range bin
    int N = 6; // position of ideal point targets
    double R[] = {10500, 11000, 12000, 12008, 13000, 13002}; // radar cross section
    double RCS[] = {1, 1, 1, 1, 1, 1};
    LFM_radar(T,B,Rmin,Rmax,R,RCS,N);
}

int main(){
    LFM_radar();
    return 0;
}
