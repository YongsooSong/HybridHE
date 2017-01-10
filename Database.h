#ifndef DATABASE_H
#define DATABASE_H

#include "params.h"
#include "distrib.h"

// return: number of distinct encoded polynomials 
int Read_data(PolyModQ**& DataPoly0, PolyModQ**& DataPoly1, char* filename);
 
#endif