#ifndef DOUBLEARRAY_H_
#define DOUBLEARRAY_H_
#include <vector>
#include <cmath>
#include <numeric>
#include "InfoOut.h"
class IntArray;
class DoubleMatrix;
/*
 * ���з�ʽ�ǰ�������
 */
class DoubleArray
{
protected:
	std::vector<double> values;
public:
	//���캯��
	DoubleArray(int n = 0) : values(n){}
	DoubleArray(const DoubleArray& src) : values(src.values){}
	DoubleArray& operator=(const DoubleArray& src) { values = src.values; return *this; }
	DoubleArray(DoubleArray&& src) : values(std::move(src.values)){}
	DoubleArray& operator=(DoubleArray&& src) { values = std::move(src.values); return *this; }
	//���õĹ��캯��
	DoubleArray(double) = delete;
	//��������
	~DoubleArray(){}
	//������ Ϊ��ѭ��
	std::vector<double>::iterator begin() { return this->values.begin(); }
	std::vector<double>::iterator end() { return this->values.end(); }
	std::vector<double>::const_iterator begin() const { return this->values.begin(); }
	std::vector<double>::const_iterator end() const { return this->values.end(); }
	//��ȡָ��͸���
	//��ȡֵ��ָ��͸���
	inline const double* GetPointer() const { return values.data(); }
	inline double* GetPointer() { return values.data(); }
	DoubleArray* GetCopy() { return new DoubleArray(*this); }
	//ȡֵ
	inline double& operator()(const int& i) { return values[i - 1]; }
	inline double operator() (const int& i) const { return values[i - 1]; }
	int getSize() const { return (int)values.size();}
	double& at(int i) { return values[i]; }
	double at(int i) const { return values[i]; }
	//�б��ʼ��
	inline DoubleArray(std::initializer_list<double> list) : values(list) {}
	inline DoubleArray& operator=(std::initializer_list<double> list) { values = list; return *this; }

	//�ж��Ƿ���ڷ�����ֵ�ͷǿ�
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
	//��ֵ
	void SetValue(const int& i, const double& val) { values[i - 1] = val; }
	void push_back(const double &val) { values.push_back(val); }
	void reserve(int s);
	//�����С
	void resize(const int s);

	//����
	void zero();

	//���
	void clear() { values.clear(); }

	// ��ֵȫ��Ϊ��ͬ��ֵ
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