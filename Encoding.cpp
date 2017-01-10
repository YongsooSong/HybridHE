
#include "Encoding.h"
#include <iostream>

using namespace std;

const uint64_t r[5][5]={ {0,36,3,41,18}, {1,44,10,45,2}, {62,6,43,15,61}, {28,55,25,21,56},{27,20,39,8,14}};

inline int mod (int a, int b){
    if(b < 0)
        return mod(-a, -b);
    int ret = a % b;
    if(ret < 0)
        ret += b;
    return ret;
}

uint64_t** sha3_round(uint64_t** A, uint64_t RC){
    int x, y;
    uint64_t C[5];
    uint64_t D[5];
    uint64_t B[5][5];

    //Theta step
    for(x = 0; x < 5; x++){C[x] = A[x][0] ^ A[x][1] ^ A[x][2] ^ A[x][3] ^ A[x][4];}
    for(x = 0; x < 5; x++){D[x] = C[(x + 4) % 5] ^ ((C[(x + 1) % 5] << 1) | (C[(x + 1) % 5] >> 63));}
    for(x = 0; x < 5; x++){for(y = 0; y < 5; y++){A[x][y] = A[x][y] ^ D[x];}}

    //Rho and phi step
    for(x = 0;x < 5; x++){for(y = 0; y < 5; y++){B[y][mod((2 * x + 3 * y),5)] = ((A[x][y] << r[x][y]) | (A[x][y] >> (64-r[x][y])));}}

    //Xi step
    for(x = 0; x < 5; x++){for(y = 0; y < 5; y++){A[x][y] = B[x][y] ^ ((~B[mod((x+1),5)][y]) & B[mod((x+2),5)][y]);}}

    //XOR step
    A[0][0] = A[0][0] ^ RC;
    return A;
}

//Round constant in keccak//
const uint64_t RC[24]={ 0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
                        0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
                        0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
                        0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
                        0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
                        0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008};

//keccak fuction with sha3 round function: input and output is 5 by 5 uint64 matrix
uint64_t **keccak_f(uint64_t **A){
  for(int i = 0; i < 24; i++){
    A = sha3_round(A, RC[i]);
  }
  return A;
}

void HASH(int64_t x, ZmodQ& deg0, ZmodQ& deg1){
    uint64_t ** A = new uint64_t*[5];
    for(int i = 0; i < 5 ; ++i){
        A[i] = new uint64_t[5];
        for(int j = 0; j < 5; ++j){
            A[i][j] = 0;
        }
    }
    A[0][0] = (uint64_t) x;
    A = keccak_f(A);
    deg0 = (ZmodQ) (A[0][0] % N);
    deg1 = (ZmodQ) (A[0][1] % N);   
}

void Encode_deg(ZmodQ& deg0, ZmodQ& deg1, char* ch, char* pos){
    int64_t ch1;
    if(ch[0] == 'X') ch1 = 0;
    else if(ch[0] == 'Y') ch1 = 23;
    else ch1 = atol(ch);

    int64_t P = 24 * atol(pos) + ch1;
    HASH(P, deg0, deg1);
}

void Encode_coef(ZmodQ*& SNP, char*& strREF, char*& strALT){
    int64_t temp, intSNP;
    int j;
    temp = 1;
    if(strREF[0] != 0x20){
        for(j = 0; j< SNPlen; ++j){
            if(strREF[j] == '\0'){ break; }
            else if(strREF[j] == 'A') temp <<= 2;
            else if(strREF[j] == 'T'){ temp <<= 2; temp += 1; }
            else if(strREF[j] == 'G'){ temp <<= 2; temp += 2; }
            else if(strREF[j] == 'C'){ temp <<= 2; temp += 3; }
        }
    }

    intSNP = temp << SNPbits;
    temp = 1;
    if(strALT[0] != 0x20){
        for(j = 0; j< SNPlen; ++j){
            if(strALT[j] == '\0'){ break; }
            else if(strALT[j] == 'A') temp <<= 2;
            else if(strALT[j] == 'T'){ temp <<= 2; temp += 1; }
            else if(strALT[j] == 'G'){ temp <<= 2; temp += 2; }
            else if(strALT[j] == 'C'){ temp <<= 2; temp += 3; }
        }
    }
    intSNP += temp;

    for(j = 0; j < nSNP; ++j){
        SNP[j] = (ZmodQ) (intSNP % T);
        intSNP >>= logT;
    }
}


 