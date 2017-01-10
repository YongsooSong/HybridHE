#include <random>
#include <iostream>
#include "Database.h"
#include "Encoding.h"

#define MAX 256

using namespace std;

int Read_data(PolyModQ**& DataPoly0, PolyModQ**& DataPoly1, char* filename){
    int i,j,k;
    char ch;
    char buf[MAX];
    int length;
    int nLine = 0;

    FILE * fp = fopen(filename, "rb");
    while(1){
        fgets(buf, MAX, fp);
        length = strlen(buf);
        if(length == 0) break;
        nLine++;
        buf[0] = '\0';
    }
    
    fseek(fp, 0L, SEEK_SET);
    
    char ** strChr = new char*[nLine];
    char ** strPOS = new char*[nLine];
    char ** strREF = new char*[nLine];
    char ** strALT = new char*[nLine];
    
    for(i = 0; i < nLine ; ++i){
        strChr[i] = new char[3];
        strPOS[i] = new char[20];
        strALT[i] = new char[SNPlen + 1];
        strREF[i] = new char[SNPlen + 1];
    }
    
    for(i = 0; i < nLine; ++i){
        //chrom
        for(j = 0; j < 3; ++j){
            ch = fgetc(fp);
            if(ch == 0x09){ strChr[i][j] = '\0'; break; }
            strChr[i][j] = ch;
        }
        if(j == 3){ cout << "Chrom Number Error" << endl; }

        //POS
        for(j = 0; j < 19; ++j){
            ch = fgetc(fp);
            if(ch == 0x09){ strPOS[i][j] = '\0'; break; }
            strPOS[i][j] = ch;
        }
        if(j == 19){ cout << "Position Error" << endl; }

        //ID
        do { ch = fgetc(fp); } while(ch != 0x09);

        //REF
        for(j = 0; j < SNPlen; ++j){
            ch = fgetc(fp);
            if(ch == 0x09){ strREF[i][j] = '\0'; break; }
            strREF[i][j] = ch;
        }
        if(j == SNPlen){
            strREF[i][SNPlen] = '\0';
            do { ch = fgetc(fp); } while(ch != 0x09);
        }

        //ALT
        for(j = 0; j < SNPlen; ++j){
            ch = fgetc(fp);
            if(ch == 0x09){ strALT[i][j] = '\0'; break;}
            strALT[i][j] = ch;
        }
        if(j == SNPlen) strALT[i][SNPlen] = '\0';

        //End Of Line
        do { ch = fgetc(fp); } while(ch != 0x0a);
    }
    fclose(fp);

    /****************************/
    // Conversion of REF & ALT
    /****************************/

//    for (i = 0; i < nLine; i++) cout << strChr[i] << '\t' << strPOS[i] << '\t' << strREF[i] << '\t' << strALT[i] << endl;

    ZmodQ** SNP = new ZmodQ*[nLine];
    
    for(i = 0; i < nLine; ++i){
        SNP[i] = new ZmodQ[nSNP];
    }

    for(i = 0; i < nLine; ++i){ Encode_coef(SNP[i], strREF[i], strALT[i]); }

    // intSNP[i] -> decompose as (logt)-bits integer


    /****************************/
    // Polynomial Encoding
    /****************************/

    bool ** check0 = new bool*[cnt];
    bool ** check1 = new bool*[cnt];
    for(i = 0; i < cnt; ++i){
        check0[i]= new bool[N];
        check1[i]= new bool[N];
        for(k = 0; k < N; ++k){
            check0[i][k] = 0;
            check1[i][k] = 0;
        }
    }

    ZmodQ* deg0 = new ZmodQ[nLine];
    ZmodQ* deg1 = new ZmodQ[nLine];

    // i=0
    Encode_deg(deg0[0], deg1[0], strChr[0], strPOS[0]);
    check0[0][deg0[0]] = 1;
    check1[0][deg1[0]] = 1;

    for(j = 0; j < nSNP; ++j){
        DataPoly0[0][j][deg0[0]] = unif_T(gen);
        DataPoly1[0][j][deg1[0]] = SNP[0][j] - DataPoly0[0][j][deg0[0]];
    }
        
    int nPoly = 1; // #(maixaml distinct polynomials)

    for(k = 1; k < nLine; ++k){
        Encode_deg(deg0[k], deg1[k], strChr[k], strPOS[k]);
        i = 0;
        while(check0[i][deg0[k]] == 1 && check1[i][deg1[k]] == 1){ i++; }
        if (i == nPoly) { nPoly++; }

        if(check0[i][deg0[k]] == 0 && check1[i][deg1[k]] == 0){
            for(j = 0; j < nSNP; ++j){
                DataPoly0[i][j][deg0[k]] = unif_T(gen);
                DataPoly1[i][j][deg1[k]] = SNP[k][j] - DataPoly0[i][j][deg0[k]];
            }
        }
        else if(check0[i][deg0[k]] == 0 && check1[i][deg1[k]] == 1){
            for(j = 0; j < nSNP; ++j){ DataPoly0[i][j][deg0[k]] = SNP[k][j] - DataPoly1[i][j][deg1[k]]; }            
        }
        else{
            for(j = 0; j < nSNP; ++j){ DataPoly1[i][j][deg1[k]] = SNP[k][j] - DataPoly0[i][j][deg0[k]]; }            
        }
        check0[i][deg0[k]] = 1; check1[i][deg1[k]] = 1;
    }

    for(i = 0; i < nPoly;++i){
        for(k = 0; k < N; ++k){
            if(check0[i][k] == 0){
                for(j = 0; j < nSNP; ++j){
                    DataPoly0[i][j][k] = unif_T(gen);
                }
            }
            if(check1[i][k] == 0){
                for(j = 0; j < nSNP;++j){
                    DataPoly1[i][j][k] = unif_T(gen);
                }
            }
        }
    }

    delete[] strChr;
    delete[] strPOS;
    delete[] strREF;
    delete[] strALT;

    delete[] SNP;

    delete[] check0;
    delete[] check1;

    delete[] deg0;
    delete[] deg1;

    cout << "nLine: " << nLine << ". nPoly: " << nPoly << ". Encrypted Database: " << (0.5*nSNP*nPoly) << " MB. Encrypted Result: " << (0.125*nSNP*nPoly) << " MB." << endl;
    return nPoly;
}

