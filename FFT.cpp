#include <fftw3.h>
#include "FFT.h"
#include "params.h"
#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace std;

double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;


void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2 * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2 * N, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2 * N, out, in,  FFTW_PATIENT);
}

void FFTforward(Ring_FFT res, const PolyModQ val)
{
    for(int k = 0; k < N; ++k){
        in[k]= (double) val[k];
        in[k+N]= 0.0;
    }
    
    fftw_execute(plan_fft_forw);
    
    for(int k = 0; k < N2; ++k){
        res[k][0]= (double)(out[2 * k + 1][0]);
        res[k][1]= (double)(out[2 * k + 1][1]);
    }  
}
 
void FFTbackward(PolyModQ res, const Ring_FFT val){
    int k2 = 0;
    
    for (int k = 0; k < N2; ++k){
        out[k2+1][0]= (double) val[k][0] / N;
        out[k2+1][1]= (double) val[k][1] / N;
        
        out[k2][0]   = (double) 0;
        out[k2][1]   = (double) 0;
        
        k2 += 2;
    }
    
    fftw_execute(plan_fft_back);
 
    for (int k = 0; k < N; ++k){ res[k] = (int64_t) round(in[k]); }
    //res[k] = (long int) round(in[k]); 
}

void coutPoly(PolyModQ& poly){
    int num = 0;
    for(int i = 0; i < N; ++i){
        if(poly[i] != 0){
            if(poly[i] >0){ cout << "+"; }
            cout << poly[i] << "X^" << i << " ";
            num++;
        }
    }
    if(num == 0){ cout << "zero poly"; }
    cout << endl;
}

void PolyAdd(PolyModQ& res, PolyModQ& poly1, PolyModQ& poly2){
    for(int i = 0; i < N; ++i){ res[i] = poly1[i] + poly2[i]; }
}

void PolyAddAndEqual(PolyModQ& res, PolyModQ& poly){
    for(int i = 0; i < N; ++i){ res[i] += poly[i]; }
}

void PolySub(PolyModQ& res, PolyModQ& poly1, PolyModQ& poly2){
    for(int i = 0; i < N; ++i){ res[i] = poly1[i] - poly2[i]; }
}    

void PolySubAndEqual(PolyModQ& res, PolyModQ& poly){
    for(int i = 0; i < N; ++i){ res[i] -= poly[i]; }
}

void ConMult(PolyModQ& poly, ZmodQ cons){
    for(int i = 0; i < N; ++i){ poly[i] *= cons; }
}

void PolyMult(PolyModQ& res, PolyModQ& poly1, PolyModQ& poly2){
    Ring_FFT Fpoly1, Fpoly2, Fpoly;
    FFTforward(Fpoly1, poly1);
    FFTforward(Fpoly2, poly2);

    for(int i = 0; i < N2; ++i){
        Fpoly[i][0]= (((double) Fpoly1[i][0]) * ((double)Fpoly2[i][0]))- (((double)Fpoly1[i][1]) * ((double)Fpoly2[i][1]));
        Fpoly[i][1]= (((double) Fpoly1[i][0]) * ((double)Fpoly2[i][1]))+ (((double)Fpoly1[i][1]) * ((double)Fpoly2[i][0]));
    }
    FFTbackward(res,Fpoly);
}

void PolyMultAndEqual(PolyModQ& res, PolyModQ& poly){
    Ring_FFT Fpoly, Fpoly1, Fpoly2;
    FFTforward(Fpoly1, res);
    FFTforward(Fpoly2, poly);

    for(int i = 0; i < N2; ++i){
        Fpoly[i][0]= (((double) Fpoly1[i][0]) * ((double)Fpoly2[i][0]))- (((double)Fpoly1[i][1]) * ((double)Fpoly2[i][1]));
        Fpoly[i][1]= (((double) Fpoly1[i][0]) * ((double)Fpoly2[i][1]))+ (((double)Fpoly1[i][1]) * ((double)Fpoly2[i][0]));
    }
    FFTbackward(res,Fpoly);
}