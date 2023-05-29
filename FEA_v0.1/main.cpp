// This code is the main code, 
// Functions include reading inp filesand performing finite element analysis

#include <iostream>
#include"FEA.h"

using namespace std;

int main()
{
	string read_path = "E:\\CADCAE\\FE_model\\Job-1.inp";
	string save_path = "E:\\CADCAE\\output\\";

	double em = 2.1e5, nu = 0.3;

	FEA fea_model;
	fea_model.Process(read_path, save_path, em, nu);

	system("pause");
	return 0;
}