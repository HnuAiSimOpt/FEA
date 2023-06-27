#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "IntMatrix.h"
#include "IntArray.h"
#include "DoubleMatrix.h"
#include "DoubleArray.h"
#include "Element.h"
#include "Matrial.h"

class PostProcessStress
{
protected:

public:
	PostProcessStress() {}  // (data: 2023 / 5 / 30)

	// calculate stress  (data: 2023 / 5 / 30)
	void GetStress(IntArray node_map, IntMatrix& connections, 
		DoubleMatrix& coordinates, DoubleArray& dis, DoubleArray& stress, ElasticMat& mat_);

};

#endif