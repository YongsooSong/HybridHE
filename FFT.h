// Code in Ducas-Micciancio Library

#ifndef FFT_H
#define FFT_H

#include "params.h"

void FFTsetup();
void FFTforward(Ring_FFT res, const PolyModQ val);
void FFTbackward(PolyModQ res, const Ring_FFT val);

void coutPoly(PolyModQ& poly);

void ConMult(PolyModQ& poly, ZmodQ cons);

void PolyAdd(PolyModQ& res, PolyModQ& poly1, PolyModQ& poly2);
void PolyAddAndEqual(PolyModQ& res, PolyModQ& poly);

void PolySub(PolyModQ& res, PolyModQ& poly1, PolyModQ& poly2);
void PolySubAndEqual(PolyModQ& res, PolyModQ& poly);

void PolyMult(PolyModQ& res, PolyModQ& poly1, PolyModQ& poly2);
void PolyMultAndEqual(PolyModQ& res, PolyModQ& poly);

#endif
