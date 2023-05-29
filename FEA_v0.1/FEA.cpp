#include"Fea.h"

// full process of FEA  (data: 2023 / 5 / 11)
void FEA::Process(string  read_path, string  save_path, double em, double nu)
{
	/*
	Reads mesh information.
	1. total number of elements
	2. total number of nodes
	3. node coordinates
	4. node connections
	5. load node set
	6. constraint node set.
	*/
	// read total number of elements and nodes
	ReadNumOfEleNode(read_path, &n_elem_, &n_node_);
	printf("The number of total elements : %d\n", n_elem_);
	printf("The number of total nodes : %d\n", n_node_);

	// read node coordinates and node connectivity
	coord_.resize(n_node_, 3);
	connect_.resize(n_elem_, 8);
	ReadInpFile(read_path, coord_, connect_);

	// read boundary and load
	ReadBoundaryLoad(read_path, con_set_, load_set_);

	/*
	set material parameters
	*/
	SetMat(em, nu);

	/*
	Assembly stiffness matrix.
	At present, only hexahedral isoparametric elements are considered.
	*/
	AssemblyStiffness();

	/*
	set the boundary conditions
	*/
	int n_dofs = global_stiff_mat_.GetRows();
	IntArray node_map = global_stiff_mat_.GetIndexFreeNode();
	BCs load;
	SetBCs(n_dofs , -100.0, load_set_, node_map, load);

	/*
	solving the system of linear equations
	*/
	DoubleArray dis(n_dofs);
	SuperLUsolver(global_stiff_mat_, load, dis);
	
	// export displacements and coordinates
	ExportDis(save_path, dis, node_map);
	ExportDis(save_path, coord_);
	ExportConnect(save_path, connect_);
	ExportData2VtkFile(save_path, 1.0e2, coord_, connect_, dis, node_map);
}

// set material  (data: 2023 / 5 / 24)
void FEA::SetMat(double em, double nu)
{
	mat_.em = em;
	mat_.nu = nu;
}

// Assembly stiffness matrix  (data: 2023 / 5 / 24)
void FEA::AssemblyStiffness()
{
	global_stiff_mat_.AssemblyStiffness(n_elem_, n_node_, mat_, con_set_, connect_, coord_);
}

// set the boundary conditions  (data: 2023 / 5 / 24)
void FEA::SetBCs(int ndofs, double scale, IntArray& load_idx, IntArray& reset_idx, BCs& load)
{
	load.SetBCs(ndofs, scale, load_idx, reset_idx);
}

