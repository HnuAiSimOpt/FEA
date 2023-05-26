#ifndef SpMatrix_h
#define SpMatrix_h

#include"DoubleArray.h"
#include"IntMatrix.h"
#include"IntArray.h"
#include <set>
#include"Element.h"

using std::vector;

class SpMatrix
{
protected:
	int nRows;
	int nColumns;
	int nz;             // total number of non-zero elements
	DoubleArray nzval;  // values of non-zero elements
	IntArray rowind;    // the row index for each element, which is equal to 'nz'
	IntArray colptr;    // record the row index of the first non-zero element in each column 
	                    // the length of colptr = 'total DOFs +1'
public:
	// ¹¹Ôìº¯Êý
	SpMatrix() : nRows(0), nColumns(0) {}
	SpMatrix(int m, int n) : nRows(m), nColumns(n) {}
	SpMatrix& operator=(const SpMatrix& src)
	{
		nRows = src.nRows;
		nColumns = src.nColumns;
		nz = src.nz;
		nzval = std::move(src.nzval);
		rowind = std::move(src.rowind);
		colptr = std::move(src.colptr);
		return *this;
	}
	~SpMatrix() {}

	void clear() { nzval.clear(); }
	int GetNZ() { return nz; }
	int GetRows() { return nRows; }
	int GetCols() { return nColumns; }

	// Create a ternary index of the global stiffness matrix  (data: 2023 / 5 / 11)
	void buildInternalStucture(int n_ele, IntMatrix& connections, IntArray& IndexFreeNode);

	// Get the DOFs set of the id-th element  (data: 2023 / 5 / 11)
	void GetEleDofs(IntArray& Edofs, int id, IntMatrix& connections, IntArray& IndexFreeNode);

	// assembly stiffness matrix  (data: 2023 / 5 / 12-15)
	void AssembleStiffMat(int n_elem, double em, double nu,
		IntMatrix& connections, IntArray& IndexFreeNode, DoubleMatrix& coordinates);

	// Extracts the colptr  (data: 2023 / 5 / 12-15)
	IntArray& GetColPtr() { return colptr; }

	// Extracts the rowind  (data: 2023 / 5 / 12-15)
	IntArray& GetRowIndex() { return rowind; }

	// Extracts the nzval  (data: 2023 / 5 / 12-15)
	DoubleArray& GetNzValues() { return nzval; }
};
#endif


