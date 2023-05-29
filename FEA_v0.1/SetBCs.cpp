#include"SetBCs.h"

void BCs::SetBCs(int ndofs, double scale, IntArray& load_idx, IntArray& reset_idx)
{
	load_vec_.resize(ndofs);
	load_vec_.zero();
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
			load_vec_(id_dof) = scale;
		}
	}
}