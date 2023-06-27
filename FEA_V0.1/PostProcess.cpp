//
#include"PostProcess.h"

void PostProcessStress::GetStress(IntArray node_map, IntMatrix& connections,
	DoubleMatrix& coordinates, DoubleArray& dis, DoubleArray& stress, ElasticMat& mat_)
{
	Elemens ele_hex;
	ele_hex.CalculateConstitutive(mat_.em, mat_.nu);
	int ne = connections.GetRow();  // number of elements
	DoubleMatrix dis_per_ele(24, 1);
	DoubleMatrix coord_per_ele(8, 3);
	int id_node;
	for (int i = 1; i < ne + 1; i++)
	{
		dis_per_ele.zero();
		for (int j = 1; j < 8 + 1; j++)
		{
			// Extraction coordinates
			coord_per_ele(j, 1) = coordinates(connections(i, j), 1);
			coord_per_ele(j, 2) = coordinates(connections(i, j), 2);
			coord_per_ele(j, 3) = coordinates(connections(i, j), 3);
			// Extraction displacement
			id_node = node_map(connections(i, j));
			if (id_node > 0)
			{
				dis_per_ele(3 * j - 2, 1) = dis(3 * id_node - 2);
				dis_per_ele(3 * j - 1, 1) = dis(3 * id_node - 1);
				dis_per_ele(3 * j, 1) = dis(3 * id_node);
			}
		}
		stress(i) = ele_hex.GetVonMisesStress(dis_per_ele, coord_per_ele);
	}
}

