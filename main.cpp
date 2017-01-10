#include <sys/time.h>
#include <iostream>
#include <string.h>

#include "params.h"
#include "Database.h"
#include "Encoding.h"
#include "FFT.h"
#include "Scheme.h"

using namespace std;

void task3(int Argc, char** Argv){
    if(Argc != 6){
        cout << "------------------------------------------------------" << endl;
        cerr << "Enter chromosome, Position, ref, alt, vcfFile\t"
        << "(e.g. $test 1 161276680 A T RCV000015246_10000.txt) \n ";
    }
    
    int i, j, k, l;

    char* QueryCh   =  Argv[1];
    char* QueryPOS  =  Argv[2];
    char* QueryREF  =  Argv[3];
    char* QueryALT  =  Argv[4];
    char* filename  =  Argv[5];
    
    struct timeval startTime, stopTime;
    double timeElapsed;

    cout << "------------------------------------------------------" << endl;
    cout << "Encoding the genome data ... " << endl;

    // REFpoly0: PolyModQ[cnt][nSNP]

    PolyModQ ** DataPoly0 = new PolyModQ*[cnt];
    PolyModQ ** DataPoly1 = new PolyModQ*[cnt];

    for(i = 0; i < cnt; ++i){
        DataPoly0[i]= new PolyModQ[nSNP];
        DataPoly1[i]= new PolyModQ[nSNP];
    }

    // return= #(distinct polynomials)
    int nPoly = Read_data(DataPoly0, DataPoly1, filename);

    //cerr << "npoly: " << nPoly<<endl;

    cout << "------------------------------------------------------" << endl;
    cerr << "Encoding the query ... " << endl;
    
	ZmodQ deg0, deg1;
    ZmodQ* QuerySNP = new ZmodQ[nSNP];

   	Encode_deg(deg0, deg1, QueryCh, QueryPOS);
 	deg0 = N - deg0;
	deg1 = N - deg1;
    Encode_coef(QuerySNP, QueryREF, QueryALT);

    cout << "------------------------------------------------------" << endl;
    cout << "FFT Setup ... " << endl;
    gettimeofday(&startTime, 0);
	FFTsetup();
    gettimeofday(&stopTime, 0);

    timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;
    timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;
	cout << "FFT Setup Time: " << timeElapsed / 1000 << " s" << endl;
	
    cout << "------------------------------------------------------" << endl;
	cout << "Key Generation ... " << endl;

    gettimeofday(&startTime, 0);
	PolyModQ sk;
	Scheme::skGen(sk);
	Scheme::RLWE_Cipher pk;
	Scheme::pkGen(pk, sk);
    gettimeofday(&stopTime, 0);

    timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;
    timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;
	cout << "Keygen Time: " << timeElapsed / 1000 << " s" << endl;

    cout << "------------------------------------------------------" << endl;
    cout << "Database (public key RLWE) Encryption ... " << endl;

	// Database Encryption (pk)
	Scheme::RLWE_Cipher** Data_CT0 = new Scheme::RLWE_Cipher*[nPoly];
	Scheme::RLWE_Cipher** Data_CT1 = new Scheme::RLWE_Cipher*[nPoly];

	for (i = 0; i < nPoly; ++i) {
	   Data_CT0[i] = new Scheme::RLWE_Cipher[nSNP];
	   Data_CT1[i] = new Scheme::RLWE_Cipher[nSNP];
	}

    gettimeofday(&startTime, 0);
	for (i = 0; i < nPoly; ++i) {
		for (j = 0; j < nSNP; ++j) {
			Scheme::RLWE_enc(Data_CT0[i][j], pk, DataPoly0[i][j]);
            Scheme::RLWE_enc(Data_CT1[i][j], pk, DataPoly1[i][j]);
		}
	}

    delete[] DataPoly0;
    delete[] DataPoly1;

	gettimeofday(&stopTime, 0);
    timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;
    timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;
    cout << "Enc Time: " << timeElapsed / 1000 << " s." << endl;
	
    cout << "------------------------------------------------------" << endl;
	cout << "Query (symmetric key GSW) Encryption ... " << endl;

	// Query Encryption (sk)
	Scheme::RLWE_Cipher* Query_CT0 = new Scheme::RLWE_Cipher[d2];
	Scheme::RLWE_Cipher* Query_CT1 = new Scheme::RLWE_Cipher[d2];

    gettimeofday(&startTime, 0);
	Scheme::GSW_enc(Query_CT0, sk, deg0);
	Scheme::GSW_enc(Query_CT1, sk, deg1);
    gettimeofday(&stopTime, 0);
    timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;
    timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;
	cout << "Query enc Time: " << timeElapsed / 1000 << endl;
    cout << "------------------------------------------------------" << endl;

    cout << "Evaluation ... " << endl;
    Scheme::RLWE_Cipher** Res_RCT0 = new Scheme::RLWE_Cipher*[nPoly];
    Scheme::RLWE_Cipher** Res_RCT1 = new Scheme::RLWE_Cipher*[nPoly];

    for(i = 0; i < nPoly; ++i){
    	Res_RCT0[i] = new Scheme::RLWE_Cipher[nSNP];
        Res_RCT1[i] = new Scheme::RLWE_Cipher[nSNP];
    }
    
    Scheme::LWE_Cipher** Res_CT = new Scheme::LWE_Cipher*[nPoly];

	for (i = 0; i < nPoly; ++i) {
		Res_CT[i] = new Scheme::LWE_Cipher[nSNP];
	}

    gettimeofday(&startTime, 0);
	for(i = 0; i < nPoly; ++i){
		for (j = 0; j < nSNP; ++j) {
			Scheme::mult(Res_RCT0[i][j], Data_CT0[i][j], Query_CT0);
            Scheme::mult(Res_RCT1[i][j], Data_CT1[i][j], Query_CT1);
			Scheme::AddAndEqual(Res_RCT0[i][j], Res_RCT1[i][j]);
			Scheme::RLWEtoLWE(Res_CT[i][j], Res_RCT0[i][j]);
    	}
	}
	
    gettimeofday(&stopTime, 0);

    delete[] Query_CT0;
    delete[] Query_CT1;

    delete[] Data_CT0;
    delete[] Data_CT1;

    delete[] Res_RCT0;
    delete[] Res_RCT1;

    timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;
	timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;
    cout << "Evaluation Time: " << timeElapsed / 1000 << " s" << endl;
    cout << "------------------------------------------------------" << endl;


    cerr << "Decryption ... " << endl;

	ZmodQ** Res = new ZmodQ*[nPoly];

    gettimeofday(&startTime, 0);
	
	for (i = 0; i < nPoly; ++i) {
		Res[i] = new ZmodQ[nSNP];
		for(j = 0; j < nSNP; ++j){ Scheme::LWE_dec(Res[i][j], sk, Res_CT[i][j]); }
	}
    gettimeofday(&stopTime, 0);
    delete[] Res_CT;

	timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;
    timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;
    cout << "Decryption Time: " << timeElapsed << " ms." << endl;

    cout << "------------------------------------------------------" << endl;

    cout << "Decoding ... " << endl;

    for(i = 0; i < nPoly; ++i){
        for(j = 0; j < nSNP; ++j){
            if(Res[i][j] != QuerySNP[j]) break;
        }
        if(j == nSNP) break;
    }
    if(i == nPoly){ cout << "no matching" << endl; }
    else{ cout << "matched" << endl; }

    delete[] Res;

    cout << "------------------------------------------------------" << endl;
}

int main(int Argc, char** Argv){
	task3(Argc, Argv);
}