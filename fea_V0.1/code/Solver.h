//
#ifndef Solver_h
#define Solver_h

#include "slu_ddefs.h"
#include"SpMatrix.h"

void SuperLUsolver(SpMatrix& A, DoubleArray& b, DoubleArray& x);

#endif