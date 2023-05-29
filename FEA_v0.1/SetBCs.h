#ifndef SETBCS_H
#define SETBCS_H

#include"IntArray.h"
#include"DoubleArray.h"

class BCs
{
public:
	DoubleArray load_vec_;
public:
	void SetBCs(int ndofs, double scale, IntArray& load_idx, IntArray& reset_idx);

	// Extracts the load_vec_  (data: 2023 / 5 / 29)
	DoubleArray& GetLoadVec() { return load_vec_; }
};

#endif