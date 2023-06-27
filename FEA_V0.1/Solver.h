#ifndef SOLVER_H
#define SOLVER_H

#ifndef SUPERLU
#define SUPERLU
#include "slu_ddefs.h"
#endif

#include "Assembly.h"
#include "SetBCs.h"

void SuperLUsolver(Assembly& A, BCs& b, DoubleArray& x);

#endif
