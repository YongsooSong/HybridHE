#include "Scheme.h"
#include <cstdlib>
#include <iostream>
#include <random>

using namespace std;

namespace Scheme {	
	void skGen(PolyModQ& sk) {
		int i, hw = 0;
		for(i = 0; i < N; ++i) sk[i] = 0;

		while(hw < 64){
			i = unif_N(gen);
			if(sk[i] == 0){
				sk[i] = (unif2(gen) << 1) - 1;
				hw++;
			}
		}
	}

	void pkGen(RLWE_Cipher& pk, PolyModQ& sk){
		for (int i = 0; i < N; ++i) { pk.a[i] = unif(gen); }
		PolyMult(pk.b, pk.a, sk);
		for(int i = 0; i < N; ++i){ pk.b[i] = Sample(Chi1) - pk.b[i]; }
	}

	void RLWE_enc(RLWE_Cipher& RLWEct, RLWE_Cipher& pk, PolyModQ& m) {
		PolyModQ v;
		for (int i = 0; i < N; ++i) { v[i] = (vi(gen) - 1) >> 1; }

		PolyMult(RLWEct.a, pk.a, v);
		PolyMult(RLWEct.b, pk.b, v);
		for(int i = 0; i < N; ++i){
			RLWEct.a[i] += Sample(Chi1);
			RLWEct.b[i] += (Sample(Chi1) + (m[i] << logu));
		}
	}

	void GSW_enc(RLWE_Cipher* GSWct, PolyModQ& sk, int deg){
		for(int i = 0; i < d2; ++i){
			for(int j=0; j < N; ++j){ GSWct[i].a[j] = unif(gen); }
			PolyMult(GSWct[i].b, GSWct[i].a, sk);
			for(int j = 0; j < N; ++j){ GSWct[i].b[j] = Sample(Chi1) - GSWct[i].b[j]; }
		}

		ZmodQ pow_B = 1;
		for(int i = 0; i < d; ++i){
			GSWct[2 * i].b[deg] -= pow_B;
			GSWct[2 * i + 1].a[deg] -= pow_B;
			pow_B <<= logB;
		}
	}

	void AddAndEqual(RLWE_Cipher& res, RLWE_Cipher& RLWEct) {
		for (int i = 0; i < N; ++i) {
			res.a[i] += RLWEct.a[i];
			res.b[i] += RLWEct.b[i];
		}
	}

	void mult(RLWE_Cipher& res, RLWE_Cipher& RLWEct, RLWE_Cipher* GSWct) {
		int i, j;
		PolyModQ a, b, poly;
		ZmodQ tmp;

		for(i = 0; i < N; ++i){
			b[i] = RLWEct.b[i];
			a[i] = RLWEct.a[i];
			res.a[i] = 0;
			res.b[i] = 0;
		}

		PolyModQ decomp_a, decomp_b;
		for(j = 0; j < d; ++j){
			for(i = 0; i < N; ++i){
				tmp = b[i] & B1;
				if(tmp >= B2){ tmp -= Bg; }
				b[i] = (b[i] - tmp) >> logB;
				decomp_b[i] = tmp;

				tmp = a[i] & B1;
				if(tmp >= B2){ tmp -= Bg; }
				a[i] = (a[i] - tmp) >> logB;
				decomp_a[i] = tmp;
			}

			PolyMult(poly, decomp_b, GSWct[2 * j].b);
			PolyAddAndEqual(res.b, poly);
			PolyMult(poly, decomp_b, GSWct[2 * j].a);
			PolyAddAndEqual(res.a, poly);

			PolyMult(poly, decomp_a, GSWct[2 * j + 1].b);
			PolyAddAndEqual(res.b, poly);
			PolyMult(poly, decomp_a, GSWct[2 * j + 1].a);
			PolyAddAndEqual(res.a, poly);			
		}
	}

	ZmodQ ModSwitch(ZmodQ a){
		ZmodQ tmp = a & MS;
		if(tmp > MS2) tmp -= MS;
		return (a - tmp);
	}

	void RLWEtoLWE(LWE_Cipher& LWEct, RLWE_Cipher& RLWEct) {
		LWEct.b = ModSwitch(RLWEct.b[0]);
		LWEct.a[0] = ModSwitch(RLWEct.a[0]);
		for (int i = 1; i < N; ++i) { LWEct.a[i] = - ModSwitch(RLWEct.a[N - i]); }
	}

	void LWE_dec(ZmodQ& m, PolyModQ& sk, LWE_Cipher& LWEct) {
		ZmodQ tmp = LWEct.b;
		for (int i = 0; i < N; ++i) { tmp += LWEct.a[i] * sk[i]; }
		m = (tmp + u2) >> logu;
		if(m < 0) m += T;
	}
}

