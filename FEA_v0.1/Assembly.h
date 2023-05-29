#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <set>
#include "DoubleMatrix.h"
#include "DoubleArray.h"
#include "IntArray.h"
#include "IntMatrix.h"
#include "Matrial.h"
#include "Element.h"

class Assembly
{
protected:
	int nRows;
	int nColumns;
	int nz;             // total number of non-zero elements
	DoubleArray nzval;  // values of non-zero elements
	IntArray rowind;    // the row index for each element, which is equal to 'nz'
	IntArray colptr;    // record the row index of the first non-zero element in each column 
						// the length of colptr = 'total DOFs +1'
	IntArray IndexFreeNode_;    // reset indes for free nodes

public:
	Assembly() {};   // (data: 2023 / 5 / 24)
	Assembly(int m, int n): nRows(m), nColumns(n) {}  // (data: 2023 / 5 / 24)

	void clear() { nzval.clear(); }       // (data: 2023 / 5 / 24)
	int GetNZ() { return nz; }            // (data: 2023 / 5 / 24)
	int GetRows() { return nRows; }       // (data: 2023 / 5 / 24)
	int GetCols() { return nColumns; }    // (data: 2023 / 5 / 24)
	void SetRowsCols(int row, int col) { nRows = row; nColumns = col; }  // (data: 2023 / 5 / 24)
	IntArray GetIndexFreeNode() { return IndexFreeNode_;  }  // (data: 2023 / 5 / 24)

	// Extracts the colptr  (data: 2023 / 5 / 29)
	IntArray& GetColPtr() { return colptr; }

	// Extracts the rowind  (data: 2023 / 5 / 29)
	IntArray& GetRowIndex() { return rowind; }

	// Extracts the nzval  (data: 2023 / 5 / 29)
	DoubleArray& GetNzValues() { return nzval; }

	// Assembly stiffness matrix  (data: 2023 / 5 / 24)
	void AssemblyStiffness(int n_ele, int n_node, ElasticMat& material,
		IntArray& constrain, IntMatrix& connections, DoubleMatrix& coordinates);

	// Set the index of the unconstrained nodes  (data: 2023 / 5 / 24)
	void SetIndexFreeNode(int n_node, IntArray& constrain);

	// Create a ternary index of the global stiffness matrix  (data: 2023 / 5 / 24)
	void BuildInternalStucture(int n_ele, IntMatrix& connections);

	// Get the DOFs set of the id-th element  (data: 2023 / 5 / 24)
	void GetEleDofs(IntArray& Edofs, int id, IntMatrix& connections);

	// Filling the value in SCS format  (data: 2023 / 5 / 24)
	void FillSparseMat(int n_ele, ElasticMat& mat_, IntMatrix& connections, DoubleMatrix& coordinates);
};

#endif
