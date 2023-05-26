#include"Fea.h"

// full process of FEA  (data: 2023 / 5 / 11)
void FEA::ProcessFea(string  save_path_dis, string save_path_coord, string save_path_connect)
{
	/*
	Reads mesh information.
	Including the total number of elements and nodes, node coordinates,
	node connections, load node set and constraint node set.
	*/
	// read total number of elements and nodes
	ReadNumOfEleNode(path_, &n_elem_, &n_node_);

	// read node coordinates and node connections
	DoubleMatrix coordinates(n_node_, 3);
	IntMatrix connections(n_elem_, 8);
	ReadInpFile(path_, coordinates, connections);

	// read boundary and load
	IntArray constrain, load;
	ReadBoundaryLoad(path_, constrain, load);
	printf("The number of total elements : %d\n", n_elem_);
	printf("The number of total nodes : %d\n", n_node_);

	/*
	Assembly stiffness matrix.
	At present, only hexahedral isoparametric elements are considered.
	*/
	// Set the index of the unconstrained nodes
	IntArray IndexFreeNode;
	SetIndexFreeNode(IndexFreeNode, constrain);

	n_fix_ = constrain.getSize();
	n_free_ = n_node_ - constrain.getSize();
	int n_dofs = 3 * n_free_;
	SpMatrix K_fea(n_dofs, n_dofs);
	Assembly(K_fea, coordinates, connections, IndexFreeNode);

	// set load
	DoubleArray load_vec(n_dofs);
	SetLoad(-100.0, load, IndexFreeNode, load_vec);
	// solver
	DoubleArray dis(n_dofs);
	SuperLUsolver(K_fea, load_vec, dis);
	
	// export displacements and coordinates
	ExportDis(save_path_dis, dis, IndexFreeNode);
	ExportDis(save_path_coord, coordinates);
	ExportConnect(save_path_connect, connections);
}


// assembly stiffness matrix  // (data: 2023 / 5 / 11)
void FEA::Assembly(SpMatrix &K_fea, DoubleMatrix& coordinates, IntMatrix& connections, IntArray& IndexFreeNode)
{
	double em = 2.1e5, nu = 0.3;  // Todo : 修改材料参数的设置方法 ！！！
	// assembly stiffness matrix, containing only unconstrained DOFs
	K_fea.buildInternalStucture(n_elem_, connections, IndexFreeNode);
	K_fea.AssembleStiffMat(n_elem_, em, nu, connections, IndexFreeNode, coordinates);
}

// Set the index of the unconstrained nodes  // (data: 2023 / 5 / 11)
void FEA::SetIndexFreeNode(IntArray& IndexFreeNode, IntArray& constrain)
{
	/*int j = 1, correct = 1;
	for (int i = 1; i < n_node_ + 1; i++)
	{
		if (!constrain.findFirstIndexOf(i))
		{
			IndexFreeNode.push_back(i - j + correct);
		}
		else
		{
			IndexFreeNode.push_back(0);
			if (j < n_fix_) { j++; }
			else { correct = 0; }
		}
	}*/
	int id_ = 1;
	for (int i = 1; i < n_node_ + 1; i++)
	{
		if (!constrain.findFirstIndexOf(i))
		{
			IndexFreeNode.push_back(id_);
			id_++;
		}
		else
		{
			IndexFreeNode.push_back(0);
		}
	}
}


// Set load  (data: 2023 / 5 / 15)
void FEA::SetLoad(double scale, IntArray& load_idx, IntArray& reset_idx, DoubleArray& load_vec)
{
	load_vec.zero();
	int n_load = load_idx.getSize();
	int idx_, reset_idx_, id_dof;
	for (int i = 1; i < n_load + 1; i++)
	{
		idx_ = load_idx(i);
		reset_idx_ = reset_idx(idx_);
		if (reset_idx_ > 0)
		{
			// id_dof = 3 * reset_idx_ - 2;  // X
			// id_dof = 3 * reset_idx_ - 1;  // Y
			// id_dof = 3 * reset_idx_ - 0;  // Z
			id_dof = 3 * reset_idx_ - 1;
			load_vec(id_dof) = scale;
		}
	}
}

// export displacement  (data: 2023 / 5 / 15)
void FEA::ExportDis(string path, DoubleArray& dis, IntArray& reset_idx)
{
	std::ofstream fout(path);
	if (!fout)
	{
		std::cerr << "fail to open file" << std::endl;
	}
	int n_node = reset_idx.getSize();
	int id = 1;
	for (int i = 1; i <= n_node; i++)
	{
		if (reset_idx(i))
		{
			fout << dis(3 * id - 2) << "    " << dis(3 * id - 1) << "    " << dis(3 * id) << std::endl;
			id++;
		}
		else
		{
			fout << 0.0 << "    " << 0.0 << "    " << 0.0 << std::endl;
		}
	}
}

// export coordinates  (data: 2023 / 5 / 15)
void FEA::ExportDis(string path, DoubleMatrix& coord)
{
	std::ofstream fout(path);
	if (!fout)
	{
		std::cerr << "fail to open file" << std::endl;
	}
	int n_node = coord.GetRow();
	int id = 1;
	for (int i = 1; i <= n_node; i++)
	{
		fout << coord(i, 1) << "    " << coord(i, 2) << "    " << coord(i, 3) << std::endl;
	}
}


// export connect  (data: 2023 / 5 / 21)
void FEA::ExportConnect(string path, IntMatrix& connect)
{
	std::ofstream fout(path);
	if (!fout)
	{
		std::cerr << "fail to open file" << std::endl;
	}
	int n_node = connect.GetRow();
	int id = 1;
	for (int i = 1; i <= n_node; i++)
	{
		fout << connect(i, 1) << "    " << connect(i, 2) << "    " << connect(i, 3) << "    " 
			<< connect(i, 4) << "    " << connect(i, 5) << "    " << connect(i, 6) << "    " 
			<< connect(i, 7) << "    " << connect(i, 8) << std::endl;
	}
}