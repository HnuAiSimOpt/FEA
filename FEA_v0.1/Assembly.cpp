#include"Assembly.h"

// Assembly stiffness matrix  (data: 2023 / 5 / 24)
void Assembly::AssemblyStiffness(int n_ele, int n_node, ElasticMat& material, 
	IntArray& constrain, IntMatrix& connections, DoubleMatrix& coordinates)
{
	// Set the index of the unconstrained nodes
	SetIndexFreeNode(n_node, constrain);
	// set shape of matrix
	int n_fix = constrain.getSize();
	int n_free = n_node - constrain.getSize();
	int n_dofs = 3 * n_free;
	SetRowsCols(n_dofs, n_dofs);
	// Create a ternary index of the global stiffness matrix
	BuildInternalStucture(n_ele, connections);
	//  Assembly stiffness matrix
	FillSparseMat(n_ele, material, connections, coordinates);
}

// Set the index of the unconstrained nodes  (data: 2023 / 5 / 24)
void Assembly::SetIndexFreeNode(int n_node, IntArray& constrain)
{
	int id_ = 1;
	for (int i = 1; i < n_node + 1; i++)
	{
		if (!constrain.findFirstIndexOf(i))
		{
			IndexFreeNode_.push_back(id_);
			id_++;
		}
		else
		{
			IndexFreeNode_.push_back(0);
		}
	}
}

// Create a ternary index of the global stiffness matrix  (data: 2023 / 5 / 24)
void Assembly::BuildInternalStucture(int n_ele, IntMatrix& connections)
{
	this->nz = 0;
	// save DOFs of all node for one element
	IntArray ele_per_dofs(24);
	int n_free_dofs = nRows;
	// saves the index of all DOFs corresponding to any one DOF
	std::vector<std::set<int>> columns(n_free_dofs);
	// add index of DOFs by connections
	for (int id = 1; id <= n_ele; id++)
	{
		// the DOFs index contained in the id-th elememnt
		GetEleDofs(ele_per_dofs, id, connections);
		// Iterate through each DOFs in ele_per_dofs
		int dof_row, dof_col;
		for (int id_dofs : ele_per_dofs)
		{
			if (id_dofs > 0)  // if id_dofs = 0, which is fixed DOF
			{
				dof_row = id_dofs - 1;
				for (int id_id_dofs : ele_per_dofs)
				{
					if (id_id_dofs > 0)
					{
						dof_col = id_id_dofs - 1;
						columns[dof_col].insert(dof_row); // Store by column (CSC)
					}
				}
			}
		}
	}
	// build rowind, colptr, and initialize nzval
	for (int id = 0; id < n_free_dofs; id++)
	{
		this->nz += columns[id].size();  // calculate number of elements
	}
	this->rowind.resize(nz);                   // rocord row index of each elements    
	this->colptr.resize(n_free_dofs + 1);      // record the row index of the first non-zero element in each column 
	int id_ = 0;
	for (int j = 0; j < n_free_dofs; j++)
	{
		colptr.at(j) = id_;               // int& at(int i) { return values[i]; }  
		for (int row : columns[j])
		{
			this->rowind.at(id_++) = row;
		}
	}
	colptr.at(n_free_dofs) = id_;  // solver usage rule
	nzval.resize(nz);
	nzval.zero();
}

// Get the DOFs set of the id-th element  (data: 2023 / 5 / 24)
void Assembly::GetEleDofs(IntArray& ele_per_dofs, int id, IntMatrix& connections)
{
	ele_per_dofs.zero();
	for (int j = 1; j < 8 + 1; j++)
	{
		int dof_id = connections(id, j);
		if (dof_id > 0)
		{
			ele_per_dofs(3 * j - 2) = 3 * IndexFreeNode_(dof_id) - 2;
			ele_per_dofs(3 * j - 1) = 3 * IndexFreeNode_(dof_id) - 1;
			ele_per_dofs(3 * j) = 3 * IndexFreeNode_(dof_id);
		}
	}
}

// Filling the value in SCS format  (data: 2023 / 5 / 24)
void Assembly::FillSparseMat(int n_ele, ElasticMat& mat_, IntMatrix& connections, DoubleMatrix& coordinates)
{
	int id_node;
	IntArray dofs_per_ele(24);
	DoubleMatrix coord_node_per_ele(8, 3);
	Elemens ele_hex;
	ele_hex.CalculateConstitutive(mat_.em, mat_.nu);
	for (int id_ele = 1; id_ele < n_ele + 1; id_ele++)
	{
		// all DOFs and coordinates of each element
		for (int id_ = 1; id_ < 8 + 1; id_++)
		{
			id_node = connections(id_ele, id_);
			dofs_per_ele(3 * id_ - 2) = 3 * IndexFreeNode_(id_node) - 2;
			dofs_per_ele(3 * id_ - 1) = 3 * IndexFreeNode_(id_node) - 1;
			dofs_per_ele(3 * id_) = 3 * IndexFreeNode_(id_node);
			coord_node_per_ele(id_, 1) = coordinates(id_node, 1);
			coord_node_per_ele(id_, 2) = coordinates(id_node, 2);
			coord_node_per_ele(id_, 3) = coordinates(id_node, 3);
		}
		// calculate element stiffness
		ele_hex.CalculateHexahedron(coord_node_per_ele);
		// assemble the element stiffness matrix into the global stiffness matrix
		int i_dof_temp, j_dof_temp;
		for (int mm = 0; mm < 24; mm++)
		{
			j_dof_temp = dofs_per_ele.at(mm);           // int& at(int i) { return values[i]; }
			if (j_dof_temp > 0)
			{
				int start = colptr.at(j_dof_temp - 1);  // row index of first non-zero at dof_temp-th column
				int t = start;
				int last_ii = this->nRows + 1;
				for (int nn = 0; nn < 24; nn++)
				{
					i_dof_temp = dofs_per_ele.at(nn);
					if (i_dof_temp > 0)
					{
						if (i_dof_temp < last_ii)
						{
							t = start;
						}
						else if (i_dof_temp > last_ii)
						{
							t++;
							ERROR("const IntArray& loc, const DoubleArray &mat", i_dof_temp);
						}
						for (; rowind.at(t) < i_dof_temp - 1; t++)  // Lower Triangular Matrix
						{
							if (t > colptr.at(j_dof_temp))          //colptr.at(j_dof_temp): row index of (dof_temp+1)-th column
							{
								ERROR("const IntArray& loc, const DoubleArray &mat", i_dof_temp);
								// In conjunction with the "else if judgment" mentioned earlier, it is a redundant judgment.
							}
						}
						this->nzval.at(t) += ele_hex.Ke.at(mm, nn);  // double& at(int i, int j) { return values[nColumns * i + j]; }
					}
				}
			}
		}
	}
}