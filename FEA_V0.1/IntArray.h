#ifndef INTARRAY_H_
#define INTARRAY_H_
#include <vector>
#include <algorithm>
#include <functional>
class IntArray
{
private:
	std::vector<int> values;
public:
	//���캯��
	IntArray(int n = 0) : values(n){}
	IntArray(const IntArray& src) : values(src.values){}
	IntArray& operator=(const IntArray& src){ values = src.values; return *this; }
	IntArray& operator=(IntArray&& src){ values = std::move(src.values); return *this; }
	IntArray(IntArray&& src) : values(std::move(src.values)) {}

	IntArray(double) = delete;
	
	//������ Ϊ��ѭ��
	std::vector<int>::iterator begin() { return this->values.begin(); }
	std::vector<int>::iterator end() { return this->values.end(); }
	std::vector<int>::const_iterator begin() const { return this->values.begin(); }
	std::vector<int>::const_iterator end() const { return this->values.end(); }
	//�б��ʼ��
	inline IntArray(std::initializer_list<int> list) : values(list) {}
	inline IntArray& operator=(std::initializer_list<int> list) { values = list; return *this; }
	
	//��������
	~IntArray(){}
	//ȡֵ
	inline int& operator()(const int& i) { return values[i - 1]; }
	inline int operator()(const int& i) const { return values[i - 1]; }
	int getSize() const { return (int)values.size(); }
	int& at(int i) { return values[i]; }
	int at(int i) const { return values[i]; }
	//��ȡֵ��ָ��͸���
	inline const int* GetPointer() const { return values.data(); }
	inline int* GetPointer() { return values.data(); }
	IntArray* GetCopy() { return new IntArray(*this); }
	// �ж��Ƿ�ǿ�
	inline bool isEmpty() const { return this->values.empty(); }
	//��ֵ
	void SetValue(const int& i, const int& val) { values[i - 1] = val; }
	void push_back(const int &val) { values.push_back(val); }
	//����
	void zero(){ std::fill(this->values.begin(), this->values.end(), 0); }
	//�ı��С
	void resize(const int &i){ this->values.resize(i); }
	// ����������С
	void reserve(int futuresize) { this->values.reserve(futuresize); }
	//���
	void clear() { values.clear(); }
	//����
	void sortmin2max(){ std::sort(values.begin(), values.end()); }
	void sortmax2min(){ std::sort(values.begin(), values.end(), std::greater<int>()); }
	//Ѱ�����ֵ����Сֵ
	int MaxValue() const { return *std::max_element(values.begin(), values.end()); }
	int MinValue() const { return *std::min_element(values.begin(), values.end()); }
	//ȥ����Ԫ�أ�����������Ԫ��
	void GetRidOfZero(IntArray& Input);
	//Ѱ�Ҷ�Ӧֵ������
	int findFirstIndexOf(int value) const
	{
		auto it = std::find(values.begin(), values.end(), value);
		if (it == values.end())
		{
			return 0;
		}
		else
		{
			return (int)(it - values.begin() + 1);
		}
	}
	bool contains(int val) const { return findFirstIndexOf(val) > 0; }
	// Ѱ����Ӧֵ������ (ֵ�洢��answer��)
	void FindIndexOf(int value, IntArray& answer) const;
	//�����С����
	void sort(); 
	void sort(IntArray& answer);
	// ����ֵ����value
	void assign(int value);
	//���
	void printfYouself() const;
};


#endif