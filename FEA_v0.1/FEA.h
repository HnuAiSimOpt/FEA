#ifndef FEA_H
#define FEA_H

#include "ReadInp.h"
#include "Assembly.h"
#include "SetBCs.h"
#include "Solver.h"
#include "WriteOutputFile.h"


class FEA
{
public:
	int n_elem_;      // total number of elements
	int n_node_;      // total number of nodes
	int n_free_;      // total number of free nodes
	int n_fix_;       // total number of fixed nodes
	
	DoubleMatrix coord_;     // coordinates of node
	IntMatrix connect_;      // connectivity of nodes
	IntArray con_set_;       // set of constraint nodes
	IntArray load_set_;      // set of load nodes
	// Elastic material parameters
	ElasticMat mat_;
	// global stiffness matrix
	Assembly global_stiff_mat_;    
public:
	FEA() :n_elem_(0), n_node_(0), n_free_(0), n_fix_(0) {}  // (data: 2023 / 5 / 24)

	// full process of FEA  (data: 2023 / 5 / 24)
	void Process(string  read_path, string  save_path, double em, double nu);

	// set material  (data: 2023 / 5 / 24)
	void SetMat(double em, double nu);

	// Assembly stiffness matrix  (data: 2023 / 5 / 24)
	void AssemblyStiffness();

	// set the boundary conditions  (data: 2023 / 5 / 24)
	void SetBCs(int ndofs, double scale, IntArray& load_idx, IntArray& reset_idx, BCs& load);


};
#endif
