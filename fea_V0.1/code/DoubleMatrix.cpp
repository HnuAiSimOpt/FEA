#include "DoubleMatrix.h"
#include "DoubleArray.h"
#include "IntArray.h"
#include "mathfem.h"
DoubleMatrix::~DoubleMatrix(){}
double DoubleMatrix::ComputeDeterminant() const
{
	if (!this->isSquare())
	{
		ERROR("Not square matrix");
	}
	if (nRows == 1)
	{
		return values[0];
	}
	else if (nRows == 2)
	{
		return (values[0] * values[3] - values[1] * values[2]);
	}
	else if (nRows == 3)
	{
		return (values[0] * values[4] * values[8] + values[3] * values[7] * values[2]
			  + values[1] * values[5] * values[6] - values[2] * values[4] * values[6]
			  - values[1] * values[3] * values[8] - values[0] * values[5] * values[7]);
	}
	else
	{
		ERROR("cannot compute the determinant of a matrix larger than 3x3");
	}
	//return 0.0;
}
DoubleMatrix DoubleMatrix::operator+(const DoubleMatrix& src) const
{
	if (nRows != src.nRows||nColumns != src.nColumns)
	{
		ERROR("cannot add the matrixs with different row or columns");
	}
	DoubleMatrix tmp(src.nRows, src.nColumns);
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			tmp(i, j) = (*this)(i, j) + src(i, j);
		}
	}
	return tmp;

}
DoubleMatrix DoubleMatrix::operator-(const DoubleMatrix &src) const
{
	if (nRows != src.nRows || nColumns != src.nColumns)
	{
		ERROR("cannot minus the matrixs with different row or columns");
	}
	DoubleMatrix tmp(src.nRows, src.nColumns);
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			tmp(i, j) = (*this)(i, j) - src(i, j);
		}
	}
	return tmp;
}

DoubleMatrix DoubleMatrix::operator*(const DoubleMatrix &src) const
{
	if (nColumns != src.nRows)
	{
		ERROR("cannot mutiply the matrixs with different row or columns");
	}
	DoubleMatrix tmp(nRows, src.nColumns);
	tmp.zero();
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			for (int k = 1; k <= nColumns; k++)
			{
				tmp(i, j) += (*this)(i, k) * src(k, j);
			}
			
		}
	}
	return tmp;
}
DoubleMatrix DoubleMatrix::operator*(const double& x) const
{
	DoubleMatrix tmp(nRows, nColumns);
	tmp.zero();
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			tmp(i, j) = (*this)(i, j)*x;
		}
	}
	return tmp;
}
void DoubleMatrix::zero()
{
	std::fill(this->values.begin(), this->values.end(), 0.0);
}
DoubleMatrix DoubleMatrix::Inverse() const
{
	DoubleMatrix tmp(nRows, nColumns);
	if (!isSquare())
	{
		ERROR("cannot inverse");
	}
	if (nRows == 2)
	{
		double det = this->ComputeDeterminant();
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		tmp(1, 1) = values[3] / det;
		tmp(1, 2) = -values[1] / det;
		tmp(2, 1) = -values[2] / det;
		tmp(2, 2) = values[0] / det;
		return tmp;
	}
	else if (nRows == 3)
	{
		
		double det = this->ComputeDeterminant();
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		tmp(1, 1) = (values[4] * values[8] - values[5] * values[7]) / det;
		tmp(1, 2) = -(values[1] * values[8] - values[7] * values[2]) / det;
		tmp(1, 3) = (values[1] * values[5] - values[4] * values[2]) / det;
		tmp(2, 1) = -(values[3] * values[8] - values[5] * values[6]) / det;
		tmp(2, 2) = (values[0] * values[8] - values[6] * values[2]) / det;
		tmp(2, 3) = -(values[0] * values[5] - values[3] * values[2]) / det;
		tmp(3, 1) = (values[3] * values[7] - values[4] * values[6]) / det;
		tmp(3, 2) = -(values[0] * values[7] - values[1] * values[6]) / det;
		tmp(3, 3) = (values[0] * values[4] - values[1] * values[3]) / det;
		return tmp;
	}
	return tmp;
}

void DoubleMatrix::beInverseOf(const DoubleMatrix& src)
{
	double det;
	if (!src.isSquare())
	{
		ERROR("cannot inverse matrix with different columns and rows");
	}
	this->resize(src.nRows, src.nColumns);
	if (nRows == 1)
	{
		(*this)(1, 1) = 1. / src(1, 1);
	}
	else if (nRows == 2)
	{
		det = src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1);
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		(*this)(1, 1) = src(2, 2) / det;
		(*this)(1, 2) = -src(1, 2) / det;
		(*this)(2, 1) = -src(2, 1)/ det;
		(*this)(2, 2) = src(1, 1) / det;
	}
	else if (nRows == 3)
	{
		det = src(1, 1) * src(2, 2) * src(3, 3) + src(1, 2) * src(2, 3) * src(3, 1) +
			src(1, 3) * src(2, 1) * src(3, 2) - src(1, 3) * src(2, 2) * src(3, 1) -
			src(2, 3) * src(3, 2) * src(1, 1) - src(3, 3) * src(1, 2) * src(2, 1);

		(*this)(1, 1) = (src(2, 2) * src(3, 3) - src(2, 3) * src(3, 2)) / det;
		(*this)(2, 1) = (src(2, 3) * src(3, 1) - src(2, 1) * src(3, 3)) / det;
		(*this)(3, 1) = (src(2, 1) * src(3, 2) - src(2, 2) * src(3, 1)) / det;
		(*this)(1, 2) = (src(1, 3) * src(3, 2) - src(1, 2) * src(3, 3)) / det;
		(*this)(2, 2) = (src(1, 1) * src(3, 3) - src(1, 3) * src(3, 1)) / det;
		(*this)(3, 2) = (src(1, 2) * src(3, 1) - src(1, 1) * src(3, 2)) / det;
		(*this)(1, 3) = (src(1, 2) * src(2, 3) - src(1, 3) * src(2, 2)) / det;
		(*this)(2, 3) = (src(1, 3) * src(2, 1) - src(1, 1) * src(2, 3)) / det;
		(*this)(3, 3) = (src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1)) / det;
	}
	else
	{
		// size >3 ... gaussian elimination - slow but safe
		//
		double piv, linkomb;
		DoubleMatrix tmp = src;
		// initialize answer to be unity matrix;
		this->zero();
		for (int i = 1; i <= nRows; i++) 
		{
			(*this)(i, i) = 1.0;
		}

		// lower triangle elimination by columns
		for (int i = 1; i < nRows; i++) 
		{

			piv = tmp(i, i);
			if (fabs(piv) < 1.e-24) 
			{
				ERROR("pivot (%d,%d) to close to small (< 1.e-24)", i, i);
			}

			for (int j = i + 1; j <= nRows; j++) 
			{
				linkomb = tmp(j, i) / tmp(i, i);
				for (int k = i; k <= nRows; k++) 
				{
					tmp(j, k) -= tmp(i, k) * linkomb;
				}

				for (int k = 1; k <= nRows; k++) {
					(*this)(j, k) -= (*this)(i, k) * linkomb;
				}
			}
		}

		// upper triangle elimination by columns
		for (int i = nRows; i > 1; i--) 
		{
			piv = tmp(i, i);
			for (int j = i - 1; j > 0; j--) 
			{
				linkomb = tmp(j, i) / piv;
				for (int k = i; k > 0; k--) 
				{
					tmp(j, k) -= tmp(i, k) * linkomb;
				}

				for (int k = nRows; k > 0; k--) 
				{
					// tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
					(*this)(j, k) -= (*this)(i, k) * linkomb;
				}
			}
		}

		// diagonal scaling
		for (int i = 1; i <= nRows; i++) 
		{
			for (int j = 1; j <= nRows; j++) {
				(*this)(i, j) /= tmp(i, i);
			}
		}
	}

}

void DoubleMatrix::beInverseAndTransposeOf(const DoubleMatrix& src)
{
	double det;
	if (!src.isSquare())
	{
		ERROR("cannot inverse matrix with different columns and rows");
	}
	this->resize(src.nRows, src.nColumns);
	if (nRows == 1)
	{
		(*this)(1, 1) = 1. / src(1, 1);
	}
	else if (nRows == 2)
	{
		det = src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1);
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		(*this)(1, 1) = src(2, 2) / det;
		(*this)(1, 2) = -src(2, 1) / det;
		(*this)(2, 1) = -src(1, 2) / det;
		(*this)(2, 2) = src(1, 1) / det;
	}
	else if (nRows == 3)
	{
		det = src(1, 1) * src(2, 2) * src(3, 3) + src(1, 2) * src(2, 3) * src(3, 1) +
			src(1, 3) * src(2, 1) * src(3, 2) - src(1, 3) * src(2, 2) * src(3, 1) -
			src(2, 3) * src(3, 2) * src(1, 1) - src(3, 3) * src(1, 2) * src(2, 1);

		(*this)(1, 1) = (src(2, 2) * src(3, 3) - src(2, 3) * src(3, 2)) / det;
		(*this)(1, 2) = (src(2, 3) * src(3, 1) - src(2, 1) * src(3, 3)) / det;
		(*this)(1, 3) = (src(2, 1) * src(3, 2) - src(2, 2) * src(3, 1)) / det;
		(*this)(2, 1) = (src(1, 3) * src(3, 2) - src(1, 2) * src(3, 3)) / det;
		(*this)(2, 2) = (src(1, 1) * src(3, 3) - src(1, 3) * src(3, 1)) / det;
		(*this)(2, 3) = (src(1, 2) * src(3, 1) - src(1, 1) * src(3, 2)) / det;
		(*this)(3, 1) = (src(1, 2) * src(2, 3) - src(1, 3) * src(2, 2)) / det;
		(*this)(3, 2) = (src(1, 3) * src(2, 1) - src(1, 1) * src(2, 3)) / det;
		(*this)(3, 3) = (src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1)) / det;
	}
	else
	{
		ERROR("exceed the size 3!");
	}
}

DoubleMatrix DoubleMatrix::Transpose() const
{
	DoubleMatrix tmp(nColumns, nRows);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			tmp(j, i) = (*this)(i, j);
		}
	}
	return tmp;
}
double DoubleMatrix::GetMatrixtrace() const
{
	if (!this->isSquare())
	{
		ERROR("Not square matrix");
	}
	double answer = 0.0;
	for (int i = 1; i <= nRows; i++)
	{
		answer += (*this)(i,i);
	}
	return answer;
}
//std::ostream& operator<<(std::ostream &out, const DoubleMatrix &src)
//{
//	for (int i = 1; i <= src.nRows; i++)
//	{
//		for (int j = 1; j <= src.nColumns; j++)
//		{
//			out << " " << src(i,j);
//		}
//		out << "\n";
//	}
//	return out;
//
//}
void DoubleMatrix::SetValues(const int i, const int j, const double value)
{
	(*this)(i, j) = value;
}
void DoubleMatrix::Add(const DoubleMatrix& src, double s)
{
	if (this->isEmpty())
	{
		this->resize(src.nRows, src.nColumns);
		this->zero();
	}
	if (nRows < src.nRows || nColumns < src.nColumns)
	{
		ERROR("cannot add the matrixs with large row or columns");
	}
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			(*this)(i, j) += src(i, j) * s;
		}
	}
}
void DoubleMatrix::minus(const DoubleMatrix& src)
{
	if (nRows != src.nRows || nColumns != src.nColumns)
	{
		ERROR("cannot minus the matrixs with different row or columns");
	}
	for (int i = 1; i <=nRows; i++)
	{
		for (int j = 1; j <=nColumns; j++)
		{
			(*this)(i, j) -= src(i, j);
		}
	}
}
void DoubleMatrix::Scale(const double& cof)
{
	for (auto& val : values)
	{
		val *= cof;
	}
}

DoubleArray DoubleMatrix::DotArray(const DoubleArray& src) const
{
	if (this->nColumns != src.getSize())
	{
		ERROR("Cannot Dot!");
	}
	DoubleArray answer(this->nRows);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			answer(i) += (*this)(i, j)*src(j);
		}
	}
	return answer;
}

void DoubleMatrix::Subtract(const DoubleMatrix& src, double scale)
{
	if (nRows != src.nRows && nColumns != src.nColumns)
	{
		ERROR("cannot subtract for mismatch size!");
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) -= src(i, j) * scale;
		}
	}
}

void DoubleMatrix::ComputeNMatrix(const DoubleArray& src, int dim)
{
	this->resize(dim, dim * src.getSize());
	for (int i = 1; i <= src.getSize(); i++)
	{
		for (int j = 1; j <= dim; j++)
		{
			(*this)(j, (i - 1) * dim + j) = src(i);
		}
	}
}

void DoubleMatrix::TProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	this->resize(a.GetColumn(), b.GetColumn());
	for (int i = 1; i <= a.GetColumn(); i++)
	{
		for (int j = 1; j <= b.GetColumn(); j++)
		{
			double tmp = 0.;
			for (int k = 1; k <= a.GetRow(); k++)
			{
				tmp += a(k, i) * b(k, j);
			}
			(*this)(i, j) = tmp;
		}
	}
	this->Scale(scale);
}

void DoubleMatrix::beProductOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.GetRow();
		this->nColumns = b.GetColumn();
		this->values.assign(a.GetRow() * b.GetColumn(), 0);
	}
	this->clear();
	this->resize(a.GetRow(), b.GetColumn());
	this->zero();
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			double sum = 0.;
			for (int k = 1; k <= a.GetColumn(); k++)
			{
				sum += a(i, k) * b(k, j);
			}
			(*this)(i, j) = sum * scale;
		}
	}
}

void DoubleMatrix::beProductTOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	DoubleMatrix tmp(b.Transpose());
	*this = a * tmp;
	this->Scale(scale);
}


void DoubleMatrix::plusProductSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.nColumns;
		this->nColumns = b.nColumns;
		this->values.assign(a.nColumns * b.nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = i; j <= nColumns; j++)
		{
			double sum = 0;
			for (int k = 1; k <= a.nRows; k++)
			{
				sum += a(k, i) * b(k, j);
			}
			(*this)(i, j) += sum * scale;
		}
	}

}



void DoubleMatrix::plusProductUnSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.nColumns;
		this->nColumns = b.nColumns;
		this->values.assign(a.nColumns * b.nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			double sum = 0;
			for (int k = 1; k <= a.nRows; k++)
			{
				sum += a(k, i) * b(k, j);
			}
			(*this)(i, j) += sum * scale;
		}
	}
}

void DoubleMatrix::symmetrized()
{
	for (int i = 2; i <= nRows; i++)
	{
		for (int j = 1; j < i; j++)
		{
			(*this)(i, j) = (*this)(j, i);
		}
	}
}

void DoubleMatrix::beUnitMatrix()
{
	if (this->nColumns != this->nRows)
	{
		ERROR("can not be unit matrix!");
	}
	this->zero();
	for (int i = 1; i <= this->nRows; i++)
	{
		(*this)(i, i) = 1.;
	}
}

void DoubleMatrix::beMatrixForm(const DoubleArray& src)
{
	
	if (src.getSize() == 9)
	{
		this->resize(3, 3);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(3, 3) = src(3);
		(*this)(2, 3) = src(4);
		(*this)(1, 3) = src(5);
		(*this)(1, 2) = src(6);
		(*this)(3, 2) = src(7);
		(*this)(3, 1) = src(8);
		(*this)(2, 1) = src(9);
	}
	else if (src.getSize() == 6)
	{
		this->resize(3, 3);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(3, 3) = src(3);
		(*this)(2, 3) = src(4);
		(*this)(1, 3) = src(5);
		(*this)(1, 2) = src(6);
		(*this)(3, 2) = src(4);
		(*this)(3, 1) = src(5);
		(*this)(2, 1) = src(6);
	}
	else if (src.getSize() == 4)
	{
		this->resize(2, 2);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(1, 2) = src(3);
		(*this)(2, 1) = src(4);
	}
	else if (src.getSize() == 3)
	{
		this->resize(2, 2);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(1, 2) = src(3);
		(*this)(2, 1) = src(3);
	}
	else
	{
		ERROR("unknown vector size!");
	}
}

void DoubleMatrix::beSubMatrix(const DoubleMatrix& src, const IntArray& a, const IntArray& b)
{
	this->resize(a.getSize(), b.getSize());
	for (int i = 1; i <= a.getSize(); i++)
	{
		for (int j = 1; j <= b.getSize(); j++)
		{
			(*this)(i, j) = src(a(i), b(j));
		}
	}
}

DoubleMatrix::DoubleMatrix(std::initializer_list<std::initializer_list<double>>mat)
{
	this->resize(mat.size(), mat.begin()->size());
	auto p = this->values.begin();
	for (auto col : mat)
	{
		for (auto x : col)
		{
			*p = x;
			p++;
		}
	}
}
DoubleMatrix& DoubleMatrix::operator=(std::initializer_list<std::initializer_list<double>> mat)
{
	this->resize(mat.size(), mat.begin()->size());
	auto p = this->values.begin();
	for (auto col : mat)
	{
		for (auto x : col)
		{
			*p = x;
			p++;
		}
	}
	return *this;
}

void DoubleMatrix::printfYourself() const
{
	printf("(%d x %d): \n", nRows, nColumns);
	if (nRows <= 250 && nColumns <= 250) {
		for (int i = 1; i <= nRows; ++i) {
			for (int j = 1; j <= nColumns && j <= 100; ++j) {
				printf("%10.3e  ", (*this)(i, j));
			}

			printf("\n");
		}
	}
	else {
		for (int i = 1; i <= nRows && i <= 20; ++i) {
			for (int j = 1; j <= nColumns && j <= 10; ++j) {
				printf("%10.3e  ", (*this)(i, j));
			}
			if (nColumns > 10) printf(" ...");
			printf("\n");
		}
		if (nRows > 20)  printf(" ...\n");
	}
}

void DoubleMatrix::copySubVectoRow(const DoubleArray& src, int r, int c)
{
	c--;
	int cols = src.getSize();
	int nr = r;
	int nc = c + cols;
	if (this->GetRow() < nr || this->GetColumn() < nc)
	{
		this->resizewithDatas(max_fem(this->GetRow(), nr), max_fem(this->GetColumn(), nc));
	}

	for (int i = 1; i <= cols; i++)
	{
		(*this)(nr, c + i) = src(i);
	}
}

void DoubleMatrix::resizewithDatas(int r, int c)
{// 旧数据会被保存的
	// 如果是原矩阵大小就不能加入
	if (r == this->nRows&&c == this->nColumns)
	{
		return;
	}

	// move 这个命令比较好用
	DoubleMatrix old(*this);
	this->nRows = r;
	this->nColumns = c;
	this->values.resize(r * c);

	int ii = min_fem(r, old.GetRow());
	int jj = min_fem(c, old.GetColumn());
	for (int i = 1; i <= ii; i++)
	{
		for (int j = 1; j <= jj; j++)
		{
			(*this)(i, j) = old(i, j);
		}
	}

}

void DoubleMatrix::rotatedwith(const DoubleMatrix& r, char mode)
{
	DoubleMatrix tmp;
	if (mode == 'n')
	{
		tmp.TProduct(r, *this);
		this->beProductOf(tmp, r);
	}
	else if (mode == 't')
	{
		tmp.beProductOf(r, *this);
		this->beProductTOf(tmp, r);
	}
	else
	{
		ERROR("unsupported mode!");
	}
}

void DoubleMatrix::setSubMatrix(const DoubleMatrix& src, int r, int c)
{
	r--;
	c--;
	int srcrow = src.GetRow();
	int srccol = src.GetColumn();
	for (int i = 1; i <= srcrow; i++)
	{
		for (int j = 1; j <= srccol; j++)
		{
			(*this)(r + i, c + j) = src(i, j);
		}
	}
}

void DoubleMatrix::assemble(const DoubleMatrix& src, const IntArray& loc)
{
	int r, c, size, ii, jj;
	r = src.GetRow();
	c = src.GetColumn();
	size = loc.getSize();
	if (r != c || r != size || c != size)
	{
		ERROR("cannot assemble! size not match");
	}
	for (int i = 1; i <= r; i++)
	{
		if ((ii = loc(i)))
		{
			for (int j = 1; j <= c; j++)
			{
				if ((jj = loc(j)))
				{
					(*this)(ii, jj) += src(i, j);
				}
			}
		}
		
	}
}

bool DoubleMatrix::solveForRhs(const DoubleArray& b, DoubleArray& answer, bool transpose)
{//@to debug
	if (!this->isSquare())
	{
		ERROR("cannot solve");
	}
	if (nRows != b.getSize())
	{
		ERROR("dimension mismatch!");
	}
	int pivRow;
	double piv, linkomb, help;
	DoubleMatrix trans;
	DoubleMatrix *K;
	if (transpose)
	{
		trans = this->Transpose();
		K = &trans;
	}
	else
	{
		K = this;
	}
	answer = b;
	for (int i = 1; i < nRows ; ++i)
	{
		piv = fabs((*K)(i, i));
		pivRow = i;
		for (int j = i + 1; j <= nRows; ++j)
		{
			if (fabs((*K)(j, i)) > piv)
			{
				pivRow = j;
				piv = fabs((*K)(j, i));
			}
		}
		if (piv < 1.0e-20)
		{
			return false;
		}
		if (pivRow != i)
		{
			for (int j = i; j <= nRows; j++)
			{
				help = (*K)(i, j);
				(*K)(i, j) = (*K)(pivRow, j);
				(*K)(pivRow, j) = help;
			}
			help = answer(i);
			answer(i) = answer(pivRow);
			answer(pivRow) = help;
		}
		for (int j = i + 1; j <= nRows; j++)
		{
			linkomb = (*K)(j, i) / (*K)(i, i);
			for (int k = i; k <= nRows; k++)
			{
				(*K)(j, k) -= (*K)(i, k) * linkomb;
			}
			answer(j) -= answer(i) * linkomb;
		}
	}
	for (int i = nRows; i >= 1; i--)
	{
		help = 0.;
		for (int j = i + 1; j <= nRows; j++)
		{
			help += (*K)(i, j) * answer(j);
		}
		answer(i) = (answer(i) - help) / (*K)(i, i);
	}
	return true;
}

bool DoubleMatrix::solveForRhs(const DoubleMatrix& b, DoubleMatrix& answer, bool transpose)
{
	/*if(!this->isSquare())
	{
		ERROR("cannot solve matrix %d * %d", nRows, nColumns);
	}
	if (nRows != b.GetRow())
	{
		ERROR("dimension mismatch");
	}*/
	int pivRow, nPs;
	double piv, linkomb, help;
	DoubleMatrix *mtrx, trans;
	if (transpose) 
	{
		trans = this->Transpose();
		mtrx = &trans;
	}
	else
	{
		mtrx = this;
	}

	nPs = b.GetColumn();
	answer = b;
	for (int i = 1; i < nRows; i++)
	{
		// find the suitable row and pivot
		piv = fabs((*mtrx)(i, i));
		pivRow = i;
		for (int j = i + 1; j <= nRows; j++)
		{
			if (fabs((*mtrx)(j, i)) > piv)
			{
				pivRow = j;
				piv = fabs((*mtrx)(j, i));
			}
		}

		if (fabs(piv) < 1.e-20)
		{
			ERROR("pivot too small, cannot solve %d by %d matrix", nRows, nColumns);
		}

		// exchange rows
		if (pivRow != i) 
		{
			for (int j = i; j <= nRows; j++) 
			{
				help = (*mtrx)(i, j);
				(*mtrx)(i, j) =(*mtrx)(pivRow, j);
				(*mtrx)(pivRow, j) = help;
			}

			for (int j = 1; j <= nPs; j++) 
			{
				help = answer(i, j);
				answer(i, j) = answer(pivRow, j);
				answer(pivRow, j) = help;
			}
		}
		if (fabs(piv) < 1.e-20) 
		{
			ERROR("cannot solve, zero pivot encountered");
		}

		for (int j = i + 1; j <= nRows; j++) 
		{
			linkomb = (*mtrx)(j, i) / (*mtrx)(i, i);
			for (int k = i; k <= nRows; k++) 
			{
				(*mtrx)(j, k) -= (*mtrx)(i, k) * linkomb;
			}

			for (int k = 1; k <= nPs; k++) 
			{
				answer(j, k) -= answer(i, k) * linkomb;
			}
		}
		
	}
	for (int i = nRows; i >= 1; i--) 
	{
		for (int k = 1; k <= nPs; k++) 
		{
			help = 0.;
			for (int j = i + 1; j <= nRows; j++) 
			{
				help += (*mtrx)(i, j) * answer(j, k);
			}

			answer(i, k) = (answer(i, k) - help) / (*mtrx)(i, i);
		}
	}
	return true;
}

void DoubleMatrix::beDiagonal(const DoubleArray& src)
{
	int n = src.getSize();
	this->resize(n, n);
	this->zero();
	for (int i = 1; i <= n; i++)
	{
		(*this)(i, i) = src(i);
	}
}

void DoubleMatrix::be_I_X_I_Matrix()
{
	this->resize(6, 6);
	this->zero();
	values[0] = values[1] = values[2] = 1.;
	values[6] = values[7] = values[8] = 1.;
	values[12] = values[13] = values[14] = 1.;
}

void DoubleMatrix::beIs()
{
	this->resize(6, 6);
	this->zero();
	values[0] = values[7] = values[14] = 1.;
	values[21] = values[28] = values[35] = 0.5;
}

void DoubleMatrix::beId()
{
	this->resize(6, 6);
	this->zero();
	values[0] = values[7] = values[14] = 2. / 3.;
	values[1] = values[2] = values[6] = values[8] = values[12] = values[13] = -1. / 3.;
	values[21] = values[28] = values[35] = 0.5;
}

void DoubleMatrix::beDyadicProductOf(const DoubleArray& a, const DoubleArray& b)
{
	this->resize(a.getSize(), b.getSize());
	for (int i = 1; i <= a.getSize(); i++)
	{
		for (int j = 1; j <= b.getSize(); j++)
		{
			(*this)(i, j) = a(i) * b(j);
		}
	}

}

void DoubleMatrix::plusDyadicSymm(const DoubleArray& a, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.getSize();
		this->nColumns = a.getSize();
		this->values.assign(this->nRows * this->nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = i; j <= nColumns; j++)
		{
			(*this)(i, j) += a(i) * a(j) * scale;
		}
	}
}

void DoubleMatrix::plusDyadicUnsymm(const DoubleArray& a, const DoubleArray& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.getSize();
		this->nColumns = b.getSize();
		this->values.assign(this->nRows * this->nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) += a(i) * b(j) * scale;
		}
	}
}

void DoubleMatrix::beSymPartOf(const DoubleMatrix& src)
{
	if (src.nRows != src.nColumns)
	{
		ERROR("src is not square!");
	}
	this->resize(src.nRows, src.nColumns);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) = 0.5 * (src(i, j) + src(j, i));
		}
	}
}

void DoubleMatrix::beAntiSymPartOf(const DoubleMatrix& src)
{
	if (src.nRows != src.nColumns)
	{
		ERROR("src is not square!");
	}
	this->resize(src.nRows, src.nColumns);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) = 0.5 * (src(i, j) - src(j, i));
		}
	}
}

void DoubleMatrix::eigProjectionOf(const DoubleMatrix& vec)
{
	if (vec.GetRow() == 2)
	{// 11 22 12
		this->resize(3, 2);
		(*this)(1, 1) = vec(1, 1) * vec(1, 1);
		(*this)(1, 2) = vec(1, 2) * vec(1, 2);

		(*this)(2, 1) = vec(2, 1) * vec(2, 1);
		(*this)(2, 2) = vec(2, 2) * vec(2, 2);

		(*this)(3, 1) = vec(1, 1) * vec(2, 1);
		(*this)(3, 2) = vec(1, 2) * vec(2, 2);
	}
	else if (vec.GetRow() == 3)
	{ // 11 22 33 23 13 12
		this->resize(6, 3);
		(*this)(1, 1) = vec(1, 1) * vec(1, 1);
		(*this)(1, 2) = vec(1, 2) * vec(1, 2);
		(*this)(1, 3) = vec(1, 3) * vec(1, 3);

		(*this)(2, 1) = vec(2, 1) * vec(2, 1);
		(*this)(2, 2) = vec(2, 2) * vec(2, 2);
		(*this)(2, 3) = vec(2, 3) * vec(2, 3);

		(*this)(3, 1) = vec(3, 1) * vec(3, 1);
		(*this)(3, 2) = vec(3, 2) * vec(3, 2);
		(*this)(3, 3) = vec(3, 3) * vec(3, 3);

		(*this)(3, 1) = vec(3, 1) * vec(3, 1);
		(*this)(3, 2) = vec(3, 2) * vec(3, 2);
		(*this)(3, 3) = vec(3, 3) * vec(3, 3);

		(*this)(4, 1) = vec(2, 1) * vec(3, 1);
		(*this)(4, 2) = vec(2, 2) * vec(3, 2);
		(*this)(4, 3) = vec(2, 3) * vec(3, 3);

		(*this)(5, 1) = vec(1, 1) * vec(3, 1);
		(*this)(5, 2) = vec(1, 2) * vec(3, 2);
		(*this)(5, 3) = vec(1, 3) * vec(3, 3);

		(*this)(6, 1) = vec(1, 1) * vec(2, 1);
		(*this)(6, 2) = vec(1, 2) * vec(2, 2);
		(*this)(6, 3) = vec(1, 3) * vec(2, 3);
	}
}


void DoubleMatrix::PlusProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.GetRow();
		this->nColumns = b.GetColumn();
		this->values.assign(a.GetRow() * b.GetColumn(), 0);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			double sum = 0.;
			for (int k = 1; k <= a.GetColumn(); k++)
			{
				sum += a(i, k) * b(k, j);
			}
			(*this)(i, j) += sum * scale;
		}
	}
}

void DoubleMatrix::copyMtrxOf(const DoubleMatrix& src, int rstart, int cstart, int rend, int cend)
{
	int nr = rend - rstart + 1;
	int nc = cend - cstart + 1;
	this->resize(nr, nc);
	this->zero();
	for (int i = 1; i <= nr; i++)
	{
		for (int j = 1; j <= nc; j++)
		{
			(*this)(i, j) = src(rstart + i - 1, cstart + j - 1);
		}
	}
}

/* void DoubleMatrix::beTanspose()
{

	for (int i = 1; i <= this->nRows; i++)
	{
		for (int j = i + 1; j <= this->nColumns; j++)
		{
			double tmp;
			tmp = (*this)(i, j);
			(*this)(i, j) = (*this)(j, i);
			(*this)(j, i) = tmp;
		}
	}
}*/