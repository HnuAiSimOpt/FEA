#ifndef DOUBLEMATRIX_H
#define DOUBLEMATRIX_H
#include <vector>
#include "InfoOut.h"
class DoubleArray;
class IntArray;
class DoubleMatrix
{
protected:
	//行数
	int nRows;
	//列数
	int nColumns;
	//值
	std::vector<double>values;
public:
	DoubleMatrix(int _nRows = 0, int _nColumns = 0) : nRows(_nRows), nColumns(_nColumns), values(_nRows*_nColumns)
	{
	}
	//拷贝构造函数
	DoubleMatrix(const DoubleMatrix& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(mat.values) {}
	DoubleMatrix(DoubleMatrix&& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(std::move(mat.values)){}
	DoubleMatrix& operator=(DoubleMatrix&& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = std::move(mat.values);
		return *this;
	}
	//赋值函数
	DoubleMatrix& operator=(const DoubleMatrix& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = mat.values;
		return *this;
	}
	~DoubleMatrix();
	//基本功能
	void SetValues(const int i, const int j, const double value);
	void SetnRows(const int i) { nRows = i; }
	void SetnColumns(const int j) { nColumns = j; }
	double GetValues(const int i, const int j) const { return (*this)(i, j); }
	int GetRow() const { return nRows; }
	int GetColumn() const { return nColumns; }
	void clear() { values.clear(); nRows = 0; nColumns = 0; }
	void resize(int rows, int columns) { nRows = rows; nColumns = columns; values.resize(nRows*nColumns); }
	//列表初始化
	DoubleMatrix(std::initializer_list<std::initializer_list<double>>src);
	DoubleMatrix& operator=(std::initializer_list<std::initializer_list<double>>src);
	////判断矩阵是否存在有
	//inline bool isFinite() const
	//{
	//	for (double val : values)
	//	{
	//		if (!std::isfinite(val))
	//		{
	//			return false;
	//		}
	//	}
	//	return true;
	//}
	//判断是否为方阵
	inline bool isSquare() const
	{
		return nRows == nColumns;
	}
	//判断是否为空
	inline bool isEmpty() const
	{
		return nRows == 0 || nColumns == 0;
	}
	//得到i行j列元素
	//可修改的左值
	inline double& operator() (const int& i, const int& j) { return values[nColumns * (i - 1) + j - 1]; }
	inline double operator() (const int& i, const int& j) const { return values[nColumns * (i - 1) + j - 1]; }
	double& at(int i, int j) { return values[nColumns * i + j]; }
	double at(int i, int j) const { return values[nColumns * i + j]; }
	//求矩阵行列式的值
	double ComputeDeterminant() const;
	//矩阵运算
	//注意不是相同维度也能相加
	void Add(const DoubleMatrix& src, double s = 1);
	void minus(const DoubleMatrix& src);
	void Scale(const double& cof);
	DoubleMatrix operator+(const DoubleMatrix &src) const;
	DoubleMatrix operator-(const DoubleMatrix &src) const;
	DoubleMatrix operator*(const DoubleMatrix &src) const;
	DoubleMatrix operator*(const double& x) const;
	DoubleMatrix Inverse() const;
	void beInverseOf(const DoubleMatrix& src);
	DoubleMatrix Transpose() const;
	double GetMatrixtrace() const;
	void zero();
	DoubleArray DotArray(const DoubleArray& src) const;
	void Subtract(const DoubleMatrix& src, double scale = 1.);

	/*a(Transpose) * b * scale*/
	void plusProductUnSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1.);

	// += a * b
	void PlusProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	//a(tansposed) * b * scale
	void TProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	// this = a * b
	void beProductOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	// this = a * b(transposed)
	void beProductTOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);

	//友元重载
	friend std::ostream& operator<<(std::ostream &out, const DoubleMatrix &src);

	//在给大小之后成为单元矩阵
	void beUnitMatrix();

	void printfYourself() const;
};
	
#endif