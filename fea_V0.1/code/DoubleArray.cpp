#include "DoubleArray.h"
#include "DoubleMatrix.h"
#include "IntArray.h"
#include "mathfem.h"


void DoubleArray::reserve(int s)
{
	this->values.reserve(s);
	this->zero();
}

void DoubleArray::resize(const int i)
{
	this->values.resize(i);
}
void DoubleArray::resizewithvalue(int k, int alloc)
{
	if (alloc > 0 && (int)this->values.capacity() < k)
	{
		this->values.reserve(k + alloc);
	}
	this->resize(k);
}

void DoubleArray::copySubVec(const DoubleArray& src, int si)
{
	si--;
	this->resizewithvalue(si + src.getSize());
	std::copy(src.begin(), src.end(), this->begin() + si);
}

void DoubleArray::zero()
{
	std::fill(this->values.begin(), this->values.end(), 0.0);
}
void DoubleArray::Add(const double &val)
{
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) += val;
	}
}
void DoubleArray::Add(const DoubleArray& src)
{
	/*if (src.isEmpty())
	{
		return;
	}*/
	if (getSize() == 0)
	{
		this->resize(src.getSize());
		this->zero();
	}
	if (this->getSize() < src.getSize())
	{
		ERROR("cannot add!");
	}
	for (int i = 1; i <= src.getSize(); i++)
	{
		(*this)(i) += src(i);
	}
}
void DoubleArray::Add(const int& i, const double &val)
{
	(*this)(i) += val;
}
void DoubleArray::Add(const double& factor, const DoubleArray& b)
{
	if (this->isEmpty())
	{
		this->resize(b.getSize());
		this->zero();
	}
	if (getSize() != b.getSize())
	{
		ERROR("cannot add");
		//return;
	}
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) += factor * b(i);
	}
}
DoubleArray DoubleArray::operator+(const DoubleArray& src)
{
	if (getSize() != src.getSize())
	{
		ERROR("cannot add");
	}
	DoubleArray tmp(src.getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		tmp(i) = (*this)(i) + src.values[i - 1];
	}
	return tmp;
}
DoubleArray DoubleArray::operator+(const double& val)
{
	DoubleArray tmp(getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		tmp(i) = (*this)(i) + val;
	}
	return tmp;
}

void DoubleArray::plus_tproduct(const DoubleMatrix& a, const DoubleArray& b, double s)
{
	if (a.GetRow() != b.getSize())
	{
		ERROR("cannot dot, size not match!");
	}
	if (this->isEmpty())
	{
		this->resize(a.GetColumn());
		this->zero();
	}
	for (auto i = 1; i <= a.GetColumn(); i++)
	{
		double tmp = 0.;
		for (auto j = 1; j <= a.GetRow(); j++)
		{
			tmp += a(j, i) * b(j);
		}
		(*this)(i) += tmp * s;
	}
}

void DoubleArray::assign(double value)
{
	for (int i = 0; i < values.size(); i++)
	{
		this->values[i] = value;
	}
}

void DoubleArray::minus(const double &val)
{
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) -= val;
	}
}
void DoubleArray::minus(const DoubleArray& src)
{
	if (src.isEmpty())
	{
		return;
	}
	if (getSize() != src.getSize())
	{
		ERROR("cannot add");
		//return;
	}
	for (int i = 1; i <= this->getSize(); i++)
	{
		(*this)(i) -= src.values[i - 1];
	}
}
void DoubleArray::minus(const int& i, const double& val)
{
	(*this)(i) -= val;
}
DoubleArray DoubleArray::operator-(const double& val)
{
	DoubleArray tmp(getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		tmp(i) = (*this)(i)- val;
	}
	return tmp;
}
DoubleArray DoubleArray::operator-(const DoubleArray& src) const
{
	if (getSize() != src.getSize())
	{
		ERROR("cannot add");
	}
	DoubleArray tmp(src.getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		tmp(i) = (*this)(i) - src.values[i - 1];
	}
	return tmp;
}
DoubleMatrix DoubleArray::operator*(const DoubleArray& src)
{
	if (getSize() != src.getSize())
	{
		ERROR("cannot dynic");
	}
	DoubleMatrix tmp(getSize(), getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		for (int j = 1; j <= getSize(); j++)
		{
			tmp.SetValues(i, j, (*this)(i) * src.values[j - 1]);
		}
	}
	return tmp;
}
double DoubleArray::Dot(const DoubleArray& src) const
{
	if (getSize() != src.getSize())
	{
		ERROR("cannot dot");
	}
	double val = 0;
	for (int i = 1; i <= getSize(); i++)
	{
		val += values[i - 1] * src.values[i - 1];
	}
	return val;
}

void DoubleArray::beMinOf(const DoubleArray& a, const DoubleArray& b)
{
	if (a.getSize() == 0)
	{
		*this = b;
		return;
	}
	if (b.getSize() == 0)
	{
		*this = a;
		return;
	}
	int n = a.getSize();
	if (n != b.getSize())
	{
		ERROR("dimension dismatch");
	}
	this->values.resize(n);
	for (int i = 1; i <= n; i++)
	{
		(*this)(i) = min_fem(a(i), b(i));
	}
}

void DoubleArray::beMaxOf(const DoubleArray& a, const DoubleArray& b)
{
	if (a.getSize() == 0)
	{
		*this = b;
		return;
	}
	if (b.getSize() == 0)
	{
		*this = a;
		return;
	}
	int n = a.getSize();
	if (n != b.getSize())
	{
		ERROR("dimension dismatch");
	}
	this->values.resize(n);
	for (int i = 1; i <= n; i++)
	{
		(*this)(i) = max_fem(a(i), b(i));
	}
}

double DoubleArray::SquareNorm()
{
	return std::inner_product(this->begin(), this->end(), this->begin(), 0.);
}
double DoubleArray::SquareNorm() const
{
	return std::inner_product(this->begin(), this->end(), this->begin(), 0.);
}

double DoubleArray::Norm()
{
	return sqrt(this->SquareNorm());
}

double DoubleArray::Normlize()
{
	double norm = sqrt(this->Dot(*this));
	if (norm < 1e-80)
	{
		ERROR("cannot Normlize");
	}
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) /= norm;
	}
	return norm;
}

void DoubleArray::Scale(double s)
{
	for (double& val : values)
	{
		val *= s;
	}
}

void DoubleArray::beScaledOf(const DoubleArray& src, double s)
{
	this->resize(src.getSize());
	for (int i = 1; i <= src.getSize(); i++)
	{
		(*this)(i) = src(i) * s;
	}
}

void DoubleArray::MtrxTDotVec(const DoubleMatrix& a, const DoubleArray& b)
{
	this->clear();
	if (a.GetRow() != b.getSize())
	{
		ERROR("dimension dismatch!");
	}
	this->resize(a.GetColumn());
	for (int j = 1; j <= a.GetColumn(); j++)
	{
		double sum = 0.;
		for (int i = 1; i <= a.GetRow(); i++)
		{
			sum += a(i, j) * b(i);
		}
		(*this)(j) = sum;
	}
}

void DoubleArray::MtrxDotVect(const DoubleMatrix& a, const DoubleArray& b)
{
	this->clear();
	if (a.GetColumn() != b.getSize())
	{
		ERROR("dimension dismatch!");
	}
	this->resize(a.GetRow());
	for (int i = 1; i <= a.GetRow(); i++)
	{
		double sum = 0.;
		for (int j = 1; j <= a.GetColumn(); j++)
		{
			sum += a(i, j) * b(j);
		}
		(*this)(i) = sum;
	}
}

double DoubleArray::Length(const DoubleArray& src) const
{
	if (values.size() != src.getSize())
	{
		ERROR("cannot compute from two vector!");
	}
	double dx = 0;
	for (int i = 1; i <= getSize(); i++)
	{
		dx += (values[i - 1] - src(i)) * (values[i - 1] - src(i));
	}

	return sqrt(dx);
}
double DoubleArray::Length(const DoubleArray* src) const
{
	return this->Length(*src);
}

double DoubleArray::LengthCoords(double x1, double x2, double y1, double y2, double z1, double z2)
{
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

void DoubleArray::beSymVectorOfStress(const DoubleMatrix& src)
{
	this->resize(6);
	(*this)(1) = src(1, 1);
	(*this)(2) = src(2, 2);
	(*this)(3) = src(3, 3);
	(*this)(4) = src(2, 3);
	(*this)(5) = src(1, 3);
	(*this)(6) = src(1, 2);
}

void DoubleArray::Subtract(const DoubleArray& a, const DoubleArray& b)
{
	if (a.getSize() != b.getSize())
	{
		ERROR("cannot subtract!");
	}
	this->resize(getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) = a(i) - b(i);
	}
}

void DoubleArray::Subtract(const DoubleArray& a)
{
	if(getSize() != a.getSize())
	{
		ERROR("size is not match! cannot subtract!");
	}
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) -= a(i);
	}
}

void DoubleArray::beDifferenceOf(const DoubleArray& a, const DoubleArray& b)
{
	if (b.getSize() != a.getSize())
	{
		ERROR("size is not match! cannot subtract!");
	}
	this->resize(a.getSize());
	for (int i = 1; i <= getSize(); i++)
	{
		(*this)(i) = a(i) - b(i);
	}
}

void DoubleArray::power(const double& e)
{
	for (double& x : *this)
	{
		x = pow(x, e);
	}
}

void DoubleArray::VectorProduct(const DoubleArray& a, const DoubleArray &v)
{
	if (a.getSize() != 3 || v.getSize() != 3)
	{
		ERROR("cannot vector product!");
	}
	this->resize(3);
	(*this)(1) = a(2) * v(3) - a(3) * v(2);
	(*this)(2) = a(3) * v(1) - a(1) * v(3);
	(*this)(3) = a(1) * v(2) - a(2) * v(1);
}

void DoubleArray::assemble(const DoubleArray& a, const IntArray& loc)
{
	if (loc.getSize() != a.getSize())
	{
		ERROR("loc and a do not match!");
	}
	for (int i = 1; i <= loc.getSize(); i++)
	{
		int ii = loc(i);
		if (ii)
		{
			(*this)(ii) += a(i);
		}
	}
}

void DoubleArray::assembleSquared(const DoubleArray& a, const IntArray& loc)
{
	if (loc.getSize() != a.getSize())
	{
		ERROR("loc and a do not match!");
	}
	for (int i = 1; i <= loc.getSize(); i++)
	{
		int ii = loc(i);
		if (ii)
		{
			(*this)(ii) += a(i) * a(i);
		}
	}
}

void DoubleArray::beSymVectorOfStrain(const DoubleMatrix& src)
{
	this->resize(6);
	(*this)(1) = src(1, 1);
	(*this)(2) = src(2, 2);
	(*this)(3) = src(3, 3);
	(*this)(4) = src(2, 3) + src(3, 2);
	(*this)(5) = src(1, 3) + src(3, 1);
	(*this)(6) = src(1, 2) + src(2, 1);
}

void DoubleArray::beVectorForm(const DoubleMatrix& src)
{
	this->resize(9);
	(*this)(1) = src(1, 1);
	(*this)(2) = src(2, 2);
	(*this)(3) = src(3, 3);
	(*this)(4) = src(2, 3);
	(*this)(5) = src(1, 3);
	(*this)(6) = src(1, 2);
	(*this)(7) = src(3, 2);
	(*this)(8) = src(3, 1);
	(*this)(9) = src(2, 1);
}

void DoubleArray::rotatewith(const DoubleMatrix& r, char mode)
{
	DoubleArray tmp;
	if (mode == 'n')
	{// r{T} * this
		tmp.MtrxTDotVec(r, *this);
	}
	else if (mode == 't')
	{// r * this
		tmp.MtrxDotVect(r, *this);
	}
	else
	{
		ERROR("unsupported rotated mode!");
	}
	*this = tmp;


}

void DoubleArray::copyRowVectorOf(const DoubleMatrix& src, int r, int cstart, int cend)
{
	if (cstart < 0 || cend > src.GetColumn() || r < 0 || r > src.GetRow() || cstart > cend)
	{
		return;
	}
	int c = cend - cstart + 1;
	this->resize(c);
	for (int i = 1; i <= c; i++)
	{
		(*this)(i) = src(r, cstart + i - 1);
	}
}

void DoubleArray::copyColumnVectorOf(const DoubleMatrix& src, int c, int rstart, int rend)
{
	if (c < 0 || c > src.GetColumn() || rstart < 0 || rend > src.GetRow() || rstart > rend)
	{
		return;
	}
	int r = rend - rstart + 1;
	this->resize(r);
	for (int i = 1; i <= r; i++)
	{
		(*this)(i) = src(rstart + i - 1, c);
	}
}

void DoubleArray::printYourself() const
{
	printf("size: %d\n", this->getSize());
	for (double x : *this)
	{
		printf("%10.3e  ", x);
	}
	printf("\n");
}

void DoubleArray::negated()
{
	for (auto i = 1; i <= values.size(); i++)
	{
		values[i - 1] = - values[i - 1];
	}
}

void DoubleArray::getFrequentcy(DoubleArray& fre)
{
	fre.resize((int)this->values.size());
	fre.zero();
	for (int i = 1; i <= this->values.size(); i++)
	{
		fre(i) = sqrt(this->values[i - 1]) / (2. * 4 * atan(1));
	}
}

bool DoubleArray::containsOnlyZeros() const
{
	for (double x : *this)
	{
		if (x != 0.0)
		{
			return false;
		}
	}
	return true;
}

void DoubleArray::printfToFile(FILE* fp)
{
	for (int i = 1; i <= this->values.size(); i++)
	{
		fprintf(fp, "%e\n",(*this)(i));
	}
}

double DoubleArray::addAllElements() const
{
	double sum = 0.0;
	for (int i = 1; i <= values.size(); i++)
	{
		sum += values[i - 1];
	}
	return sum;
}
