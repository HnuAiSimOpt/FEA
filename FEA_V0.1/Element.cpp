
#include"Element.h"

// build constitutive matrix  (data: 2023 / 5 / 12)
void Elemens::CalculateConstitutive(double em, double nu)
{
	Ce.resize(6, 6); // https://zhuanlan.zhihu.com/p/97700824
	Ce.zero();
	double factor = (em * (1.0 - nu)) / ((1.0 + nu) * (1.0 - 2.0 * nu));
	double a = nu / (1.0 - nu);
	double b = (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
	Ce(1, 1) = factor * 1.0;  Ce(1, 2) = factor * a;    Ce(1, 3) = factor * a;
	Ce(2, 1) = factor * a;    Ce(2, 2) = factor * 1.0;  Ce(2, 3) = factor * a;
	Ce(3, 1) = factor * a;    Ce(3, 2) = factor * a;    Ce(3, 3) = factor * 1.0;
	Ce(4, 4) = factor * b;
	Ce(5, 5) = factor * b;
	Ce(6, 6) = factor * b;
}

// build element stiffness matrix  (data: 2023 / 5 / 12-14)
void Elemens::CalculateHexahedron(DoubleMatrix& node_coord)
{
	// employing two - point Gaussian integral
	double gp_values = 1. / sqrt(3.);
	DoubleMatrix gps(8, 3);
	gps(1, 1) = -gp_values;  gps(1, 2) = -gp_values;  gps(1, 3) = -gp_values;
	gps(2, 1) =  gp_values;  gps(2, 2) = -gp_values;  gps(2, 3) = -gp_values;
	gps(3, 1) =  gp_values;  gps(3, 2) =  gp_values;  gps(3, 3) = -gp_values;
	gps(4, 1) = -gp_values;  gps(4, 2) =  gp_values;  gps(4, 3) = -gp_values;
	gps(5, 1) = -gp_values;  gps(5, 2) = -gp_values;  gps(5, 3) =  gp_values;
	gps(6, 1) =  gp_values;  gps(6, 2) = -gp_values;  gps(6, 3) =  gp_values;
	gps(7, 1) =  gp_values;  gps(7, 2) =  gp_values;  gps(7, 3) =  gp_values;
	gps(8, 1) = -gp_values;  gps(8, 2) =  gp_values;  gps(8, 3) =  gp_values;
	DoubleMatrix dN_drst(3, 8), dN_dxyz(3, 8), strain_mat(6, 24), temp_1, temp_2;
	DoubleMatrix jacobi, inv_jacobi;
	// N = (1+s0*s)(1+t0*t)
	// (s0, t0) = [[-1, -1], [1, -1], [1, 1], [-1, 1]];
	double r, s, t;
	double det_jacobi;
	Ke.resize(24, 24);
	Ke.zero();
	for (int i = 0; i < 8; i++)
	{
		r = gps(i + 1, 1); s = gps(i + 1, 2); t = gps(i + 1, 3);

		dN_drst(1, 1) = -0.125 * (1.0 - s) * (1.0 - t);
		dN_drst(1, 2) =  0.125 * (1.0 - s) * (1.0 - t);
		dN_drst(1, 3) =  0.125 * (1.0 + s) * (1.0 - t);
		dN_drst(1, 4) = -0.125 * (1.0 + s) * (1.0 - t);
		dN_drst(1, 5) = -0.125 * (1.0 - s) * (1.0 + t);
		dN_drst(1, 6) =  0.125 * (1.0 - s) * (1.0 + t);
		dN_drst(1, 7) =  0.125 * (1.0 + s) * (1.0 + t);
		dN_drst(1, 8) = -0.125 * (1.0 + s) * (1.0 + t);

		dN_drst(2, 1) = -0.125 * (1.0 - r) * (1.0 - t);
		dN_drst(2, 2) = -0.125 * (1.0 + r) * (1.0 - t);
		dN_drst(2, 3) =  0.125 * (1.0 + r) * (1.0 - t);
		dN_drst(2, 4) =  0.125 * (1.0 - r) * (1.0 - t);
		dN_drst(2, 5) = -0.125 * (1.0 - r) * (1.0 + t);
		dN_drst(2, 6) = -0.125 * (1.0 + r) * (1.0 + t);
		dN_drst(2, 7) =  0.125 * (1.0 + r) * (1.0 + t);
		dN_drst(2, 8) =  0.125 * (1.0 - r) * (1.0 + t);

		dN_drst(3, 1) = -0.125 * (1.0 - r) * (1.0 - s);
		dN_drst(3, 2) = -0.125 * (1.0 + r) * (1.0 - s);
		dN_drst(3, 3) = -0.125 * (1.0 + r) * (1.0 + s);
		dN_drst(3, 4) = -0.125 * (1.0 - r) * (1.0 + s);
		dN_drst(3, 5) =  0.125 * (1.0 - r) * (1.0 - s);
		dN_drst(3, 6) =  0.125 * (1.0 + r) * (1.0 - s);
		dN_drst(3, 7) =  0.125 * (1.0 + r) * (1.0 + s);
		dN_drst(3, 8) =  0.125 * (1.0 - r) * (1.0 + s);

		jacobi.beProductOf(dN_drst, node_coord);    // 3x8
		// strain matrix
		inv_jacobi.beInverseOf(jacobi);             // 3x3
		dN_dxyz.beProductOf(inv_jacobi, dN_drst);   // 3x8
		strain_mat.zero();
		for (int i = 0; i < 8; i++)
		{
			strain_mat(1, 3 * i + 1) = dN_dxyz(1, i + 1);
			strain_mat(2, 3 * i + 2) = dN_dxyz(2, i + 1);
			strain_mat(3, 3 * i + 3) = dN_dxyz(3, i + 1);
			strain_mat(4, 3 * i + 1) = dN_dxyz(2, i + 1);  strain_mat(4, 3 * i + 2) = dN_dxyz(1, i + 1);
			strain_mat(5, 3 * i + 2) = dN_dxyz(3, i + 1);  strain_mat(5, 3 * i + 3) = dN_dxyz(2, i + 1);
			strain_mat(6, 3 * i + 1) = dN_dxyz(3, i + 1);  strain_mat(6, 3 * i + 3) = dN_dxyz(1, i + 1);
		}
		det_jacobi = jacobi(1, 1) * jacobi(2, 2) * jacobi(3, 3) 
			+ jacobi(1, 2) * jacobi(2, 3) * jacobi(3, 1) 
			+ jacobi(1, 3) * jacobi(2, 1) * jacobi(3, 2) 
			- jacobi(1, 3) * jacobi(2, 2) * jacobi(3, 1) 
			- jacobi(2, 3) * jacobi(3, 2) * jacobi(1, 1) 
			- jacobi(3, 3) * jacobi(1, 2) * jacobi(2, 1);
		temp_1.TProduct(strain_mat, Ce);
		temp_2.beProductOf(temp_1, strain_mat);
		Ke.Add(temp_2, det_jacobi);
	}
}


// calulate von Mises stress  (data: 2023 / 5 / 30)
double Elemens::GetVonMisesStress(DoubleMatrix& dis, DoubleMatrix& coord)
{
	double r = 0, s = 0, t = 0;
	double det_jacobi;
	DoubleMatrix dN_drst(3, 8), dN_dxyz(3, 8), strain_mat(6, 24), stress, temp;
	DoubleMatrix jacobi, inv_jacobi, dN;
	//
	dN_drst(1, 1) = -0.125 * (1.0 - s) * (1.0 - t);
	dN_drst(1, 2) = 0.125 * (1.0 - s) * (1.0 - t);
	dN_drst(1, 3) = 0.125 * (1.0 + s) * (1.0 - t);
	dN_drst(1, 4) = -0.125 * (1.0 + s) * (1.0 - t);
	dN_drst(1, 5) = -0.125 * (1.0 - s) * (1.0 + t);
	dN_drst(1, 6) = 0.125 * (1.0 - s) * (1.0 + t);
	dN_drst(1, 7) = 0.125 * (1.0 + s) * (1.0 + t);
	dN_drst(1, 8) = -0.125 * (1.0 + s) * (1.0 + t);
	//
	dN_drst(2, 1) = -0.125 * (1.0 - r) * (1.0 - t);
	dN_drst(2, 2) = -0.125 * (1.0 + r) * (1.0 - t);
	dN_drst(2, 3) = 0.125 * (1.0 + r) * (1.0 - t);
	dN_drst(2, 4) = 0.125 * (1.0 - r) * (1.0 - t);
	dN_drst(2, 5) = -0.125 * (1.0 - r) * (1.0 + t);
	dN_drst(2, 6) = -0.125 * (1.0 + r) * (1.0 + t);
	dN_drst(2, 7) = 0.125 * (1.0 + r) * (1.0 + t);
	dN_drst(2, 8) = 0.125 * (1.0 - r) * (1.0 + t);
	//
	dN_drst(3, 1) = -0.125 * (1.0 - r) * (1.0 - s);
	dN_drst(3, 2) = -0.125 * (1.0 + r) * (1.0 - s);
	dN_drst(3, 3) = -0.125 * (1.0 + r) * (1.0 + s);
	dN_drst(3, 4) = -0.125 * (1.0 - r) * (1.0 + s);
	dN_drst(3, 5) = 0.125 * (1.0 - r) * (1.0 - s);
	dN_drst(3, 6) = 0.125 * (1.0 + r) * (1.0 - s);
	dN_drst(3, 7) = 0.125 * (1.0 + r) * (1.0 + s);
	dN_drst(3, 8) = 0.125 * (1.0 - r) * (1.0 + s);
	// jacobi matrix
	jacobi.beProductOf(dN_drst, coord);         // 3x8
	// strain matrix
	inv_jacobi.beInverseOf(jacobi);             // 3x3
	dN_dxyz.beProductOf(inv_jacobi, dN_drst);   // 3x8
	strain_mat.zero();
	for (int i = 0; i < 8; i++)
	{
		strain_mat(1, 3 * i + 1) = dN_dxyz(1, i + 1);
		strain_mat(2, 3 * i + 2) = dN_dxyz(2, i + 1);
		strain_mat(3, 3 * i + 3) = dN_dxyz(3, i + 1);
		strain_mat(4, 3 * i + 1) = dN_dxyz(2, i + 1);  strain_mat(4, 3 * i + 2) = dN_dxyz(1, i + 1);
		strain_mat(5, 3 * i + 2) = dN_dxyz(3, i + 1);  strain_mat(5, 3 * i + 3) = dN_dxyz(2, i + 1);
		strain_mat(6, 3 * i + 1) = dN_dxyz(3, i + 1);  strain_mat(6, 3 * i + 3) = dN_dxyz(1, i + 1);
	}
	temp.beProductOf(strain_mat, dis);   // 6x1
	stress.beProductOf(Ce, temp);        // 6x1
	// von Mises stress
	double von_Mises_stress = 0.5 * (
		(stress(1, 1) - stress(2, 1)) * (stress(1, 1) - stress(2, 1)) +
		(stress(2, 1) - stress(3, 1)) * (stress(2, 1) - stress(3, 1)) + 
		(stress(1, 1) - stress(3, 1)) * (stress(1, 1) - stress(3, 1)) + 
		6 * (stress(4, 1) * stress(4, 1) + stress(5, 1) * stress(5, 1) + stress(6, 1) * stress(6, 1))
		);
	von_Mises_stress = sqrt(von_Mises_stress);
	return von_Mises_stress;
}

// print self  (data: 2023 / 5 / 14)
void Elemens::show()
{
	int row = Ke.GetRow();
	int col = Ke.GetColumn();
	double value;
	for (int i = 1; i < row + 1; i++)
	{
		for (int j = 1; j < col + 1; j++)
		{
			value = Ke(i, j);
			printf("(%d, %d) : %f\n", i, j, value);
		}
	}
}