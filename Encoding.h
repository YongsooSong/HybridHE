#ifndef Encoding
#define Encoding

#include "params.h"

//! @ Input: ch & pos
//! @ Output: deg0, deg1
void Encode_deg(ZmodQ& deg0, ZmodQ& deg1, char* ch, char* pos);

//! @ Input: S
//! @ Output: length of output SNP
void Encode_coef(ZmodQ*& SNP, char*& strREF, char*& strALT);
 
#endif