#ifndef IntMatrix_H_
#define IntMatrix_H_

#include<vector>
class IntMatrix
{
protected:
	//����
	int nRows;
	//����
	int nColumns;
	//ֵ
	std::vector<int>values;
public:
	IntMatrix(int _nRows = 0, int _nColumns = 0) : nRows(_nRows), nColumns(_nColumns), values(_nRows*_nColumns)
	{
	}
	//�������캯��
	IntMatrix(const IntMatrix& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(mat.values) {}
	//��ֵ����
	IntMatrix& operator=(const IntMatrix& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = mat.values;
		return *this;
	}
	~IntMatrix(){}
	//�õ�i��j��Ԫ��
	//���޸ĵ���ֵ
	inline int& operator() (const int& i, const int& j) { return values[nColumns * (i - 1) + j - 1]; }
	inline int operator() (const int& i, const int& j) const { return values[nColumns * (i - 1) + j - 1]; }
	//��������
	void SetValues(const int i, const int j, const int value) { (*this)(i, j) = value; }
	void SetnRows(const int i) { nRows = i; }
	void SetnColumns(const int j) { nColumns = j; }
	int GetValues(const int i, const int j) const { return (*this)(i, j); }
	int GetRow() const { return nRows; }
	int GetColumn() const { return nColumns; }
	void clear() { values.clear(); nRows = 0; nColumns = 0; }
	void resize(int rows, int columns) { nRows = rows; nColumns = columns; values.resize(nRows*nColumns); }
	void zero() { std::fill(this->values.begin(), this->values.end(), 0); }
	
};
#endif