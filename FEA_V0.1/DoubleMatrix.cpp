#include "DoubleMatrix.h"
#include "DoubleArray.h"
#include "IntArray.h"

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