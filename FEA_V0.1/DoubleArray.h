#ifndef DOUBLEARRAY_H_
#define DOUBLEARRAY_H_
#include <vector>
#include <cmath>
#include <numeric>
#include "InfoOut.h"
class IntArray;
class DoubleMatrix;
/*
 * 排列方式是按行排列
 */
class DoubleArray
{
protected:
	std::vector<double> values;
public:
	//构造函数
	DoubleArray(int n = 0) : values(n){}
	DoubleArray(const DoubleArray& src) : values(src.values){}
	DoubleArray& operator=(const DoubleArray& src) { values = src.values; return *this; }
	DoubleArray(DoubleArray&& src) : values(std::move(src.values)){}
	DoubleArray& operator=(DoubleArray&& src) { values = std::move(src.values); return *this; }
	//禁用的构造函数
	DoubleArray(double) = delete;
	//析构函数
	~DoubleArray(){}
	//迭代器 为了循环
	std::vector<double>::iterator begin() { return this->values.begin(); }
	std::vector<double>::iterator end() { return this->values.end(); }
	std::vector<double>::const_iterator begin() const { return this->values.begin(); }
	std::vector<double>::const_iterator end() const { return this->values.end(); }
	//获取指针和复制
	//获取值的指针和复制
	inline const double* GetPointer() const { return values.data(); }
	inline double* GetPointer() { return values.data(); }
	DoubleArray* GetCopy() { return new DoubleArray(*this); }
	//取值
	inline double& operator()(const int& i) { return values[i - 1]; }
	inline double operator() (const int& i) const { return values[i - 1]; }
	int getSize() const { return (int)values.size();}
	double& at(int i) { return values[i]; }
	double at(int i) const { return values[i]; }
	//列表初始化
	inline DoubleArray(std::initializer_list<double> list) : values(list) {}
	inline DoubleArray& operator=(std::initializer_list<double> list) { values = list; return *this; }

	//判断是否存在非有限值和非空
	inline bool isFinite() const
	{
		for (double val : values)
		{
			if (!std::isfinite(val))
			{
				return false;
			}
		}
		return true;
	}
	inline bool isEmpty() const { return this->values.empty(); }
	//赋值
	void SetValue(const int& i, const double& val) { values[i - 1] = val; }
	void push_back(const double &val) { values.push_back(val); }
	void reserve(int s);
	//重设大小
	void resize(const int s);

	//置零
	void zero();

	//清空
	void clear() { values.clear(); }

	// 赋值全部为相同的值
	void assign(double value);

	// Norm
	double SquareNorm();
	double SquareNorm() const;

	// scaled
	void Scale(double s);
	// this = src * s
	void beScaledOf(const DoubleArray& src, double s);

	// output
	void printYourself() const;
};


#endif