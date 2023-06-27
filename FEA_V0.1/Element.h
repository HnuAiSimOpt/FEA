#ifndef ELEMENT_H
#define ELEMENT_H

#include"DoubleMatrix.h"

class Elemens
{
protected:
	DoubleMatrix Ce;  // constitutive matrix
	char type;
public:
	DoubleMatrix Ke;  // element stiffness matrix

	// build constitutive matrix  (data: 2023 / 5 / 12)
	void CalculateConstitutive(double em, double nu);

	// build element stiffness matrix  (data: 2023 / 5 / 12-14)
	void CalculateHexahedron(DoubleMatrix& node_coord);

	// calulate von Mises stress  (data: 2023 / 5 / 30)
	double GetVonMisesStress(DoubleMatrix& dis, DoubleMatrix& coord);

	// print self  (data: 2023 / 5 / 14)
	void show();
};

#endif
