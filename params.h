#ifndef PARAM_H
#define PARAM_H

#include <complex.h>
#include <stdint.h>
#include <fftw3.h>

typedef int32_t ZmodQ;
 
//! HE parameters
#define logN 11
#define N 2048
//! Parameters for FFT
#define N2 1025

#define logT 11
#define T 2048 // plaintext space

//! encryption constant
#define logu 21 // 32-logT
#define u 0x0200000
#define u2 0x0100000
#define MS 0xffff
#define MS2 0x7fff

#define cnt 100  // maximal number of polynomials: 3 for 10K, 5 for 100K

//! @ DNA parameters
#if 1
#define SNPlen 2	 	// length of SNP that we will consider
#define SNPbits 5	 	// 2*SNPlen + 1 // bits of SNPs
#define nSNP 1 		 	// SNPlen2/logT // REF, ALT - > nSNP partition

#else
#define SNPlen 5	 	// length of SNP that we will consider
#define SNPbits 11	 	// 2*SNPlen + 1 // bits of SNPs
#define nSNP 2 		 	// SNPlen2/logT // REF, ALT - > nSNP partition
#endif

//! decomposition bound
#define logB 7
#define Bg 128
#define B1 127
#define B2 64

//! decomposition exponent
#define d 5 // (31/logB) + 1;
#define d2 10

//! polynomial in ZmodQ[x]/(x^n+1)
typedef ZmodQ PolyModQ[N];

typedef fftw_complex Ring_FFT[N2];
#endif
