#ifndef FEA_H
#define FEA_H

#include "ReadInp.h"
#include "SpMatrix.h"
#include"Solver.h"

class FEA
{
public:
	string path_;
	int n_elem_;     // total number of elements
	int n_node_;     // total number of nodes
	int n_free_;     // total number of free nodes
	int n_fix_;      // total number of fixed nodes
	double em_;      // elastic modulus
	double nu_;      // Poisson's ratio
public:
	FEA() :n_elem_(0), n_node_(0), n_fix_(0){}  // (data: 2023 / 5 / 11)
	FEA(string path) :n_elem_(0), n_node_(0), n_fix_(0), path_(path){}  // (data: 2023 / 5 / 11)

	// full process of FEA  (data: 2023 / 5 / 11)
	void ProcessFea(string  save_path_dis, string save_path_coord, string save_path_connect);

	// assembly stiffness matrix  // (data: 2023 / 5 / 11)
	void Assembly(SpMatrix& K_fea, DoubleMatrix& coordinates, IntMatrix& connections, IntArray& IndexFreeNode);

	// set matrial parameters
	void SetMatrial(double em, double nu) 
	{
		em_ = em;
		nu_ = nu;
	}

	// Set the index of the unconstrained nodes  (data: 2023 / 5 / 11)
	void SetIndexFreeNode(IntArray& IndexFreeNode, IntArray& constrain);


	// Set load  (data: 2023 / 5 / 15)
	void SetLoad(double scale, IntArray& load_idx, IntArray& reset_idx, DoubleArray& load_vec);

	// export displacement  (data: 2023 / 5 / 15)
	void ExportDis(string path, DoubleArray& dis, IntArray& reset_idx);

	// export coordinates  (data: 2023 / 5 / 15)
	void ExportDis(string path, DoubleMatrix& coord);

	// export connect  (data: 2023 / 5 / 21)
	void ExportConnect(string path, IntMatrix& connect);

};
#endif