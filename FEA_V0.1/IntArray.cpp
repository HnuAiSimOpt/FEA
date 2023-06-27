#include "IntArray.h"

void IntArray::sort()
{
	std::sort(this->begin(), this->end());
}
void IntArray::sort(IntArray& answer)
{
	answer.clear();
	answer = *this;
	answer.sort();
}

void IntArray::assign(int value)
{
	for (int i = 1; i <= values.size(); i++)
	{
		(*this)(i) = value;
	}
}

void IntArray::printfYouself() const
{
	printf("size: %d\n", this->getSize());
	for (int x : *this)
	{
		printf("%5d", x);
	}
	printf("\n");
}

void IntArray::GetRidOfZero(IntArray& Input)
{
	int size = 0;
	for (const int& x : Input.values)
	{
		if (x)
		{
			size++;
		}
	}
	int pos = 1;
	this->values.resize(size);
	for (int i = 1; i <= Input.getSize(); i++)
	{
		if (Input(i))
		{
			(*this)(pos++) = i;
		}
	}
}

void IntArray::FindIndexOf(int value, IntArray& answer) const
{
	for (int i = 1; i <= values.size(); i++)
	{
		if (value == values[i - 1])
		{
			answer.push_back(i);
		}
	}
}
