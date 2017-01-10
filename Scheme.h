#ifndef SCHEME_H_
#define SCHEME_H_

#include "params.h"
#include "distrib.h"
#include "FFT.h"

using namespace std;

namespace Scheme {

	typedef struct { // n-dimensional LWE ciphertext modulo q
		ZmodQ a[N]; // array of n integers modulo q 
		ZmodQ b;    // integer modulo q
	} LWE_Cipher;

	typedef struct { // n-dimensional LWE ciphertext modulo q
		PolyModQ a; // array of n integers modulo q 
		PolyModQ b;    // integer modulo q
	} RLWE_Cipher;

  	void skGen(PolyModQ& sk);
  	void pkGen(RLWE_Cipher& pk, PolyModQ& sk);
	
	void RLWE_enc(RLWE_Cipher& RLWEct, RLWE_Cipher& pk, PolyModQ& m);
	void AddAndEqual(RLWE_Cipher& res, RLWE_Cipher& RLWEct);
	
	void GSW_enc(RLWE_Cipher* GSWct, PolyModQ& sk, int deg);
	void mult(RLWE_Cipher& res, RLWE_Cipher& RLWEct, RLWE_Cipher* GSWct);
	
	void RLWEtoLWE(LWE_Cipher& LWEct, RLWE_Cipher& RLWEct);
	void LWE_dec(ZmodQ& m, PolyModQ& sk, LWE_Cipher& LWEct);

	ZmodQ ModSwitch(ZmodQ a);
};

#endif /* SCHEME_H_ */
