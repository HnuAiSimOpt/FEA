#include "DoubleArray.h"
#include "DoubleMatrix.h"
#include "IntArray.h"

void DoubleArray::reserve(int s)
{
	this->values.reserve(s);
	this->zero();
}

void DoubleArray::resize(const int i)
{
	this->values.resize(i);
}

void DoubleArray::zero()
{
	std::fill(this->values.begin(), this->values.end(), 0.0);
}

void DoubleArray::assign(double value)
{
	for (int i = 0; i < values.size(); i++)
	{
		this->values[i] = value;
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

void DoubleArray::printYourself() const
{
	printf("size: %d\n", this->getSize());
	for (double x : *this)
	{
		printf("%10.3e  ", x);
	}
	printf("\n");
}

