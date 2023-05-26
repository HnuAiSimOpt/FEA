// This code is the main code, 
// Functions include reading inp filesand performing finite element analysis

#include <iostream>
#include"Fea.h"

using namespace std;

int main()
{
	//Elemens ele_hex;
	//double em = 1.0, nu = 0.3;
	//ele_hex.CalculateConstitutive(em, nu);
	//DoubleMatrix coord_node_per_ele(8, 3);
	//coord_node_per_ele(1, 1) = 0;   coord_node_per_ele(1, 2) = 0;   coord_node_per_ele(1, 3) = 0;
	//coord_node_per_ele(2, 1) = 1;   coord_node_per_ele(2, 2) = 0;   coord_node_per_ele(2, 3) = 0;
	//coord_node_per_ele(3, 1) = 1;   coord_node_per_ele(3, 2) = 1;   coord_node_per_ele(3, 3) = 0;
	//coord_node_per_ele(4, 1) = 0;   coord_node_per_ele(4, 2) = 1;   coord_node_per_ele(4, 3) = 0;
	//coord_node_per_ele(5, 1) = 0;   coord_node_per_ele(5, 2) = 0;   coord_node_per_ele(5, 3) = 1;
	//coord_node_per_ele(6, 1) = 1;   coord_node_per_ele(6, 2) = 0;   coord_node_per_ele(6, 3) = 1;
	//coord_node_per_ele(7, 1) = 1;   coord_node_per_ele(7, 2) = 1;   coord_node_per_ele(7, 3) = 1;
	//coord_node_per_ele(8, 1) = 0;   coord_node_per_ele(8, 2) = 1;   coord_node_per_ele(8, 3) = 1;
	//ele_hex.CalculateHexahedron(coord_node_per_ele);
	//ele_hex.show();


	string path = "E:\\CADCAE\\FE_model\\Job-1.inp";
	string save_path_dis = "E:\\CADCAE\\dis.txt";
	string save_path_coord = "E:\\CADCAE\\coord.txt";
	string save_path_connect = "E:\\CADCAE\\connect.txt";
	FEA fea_model(path);
	fea_model.ProcessFea(save_path_dis, save_path_coord, save_path_connect);

	double em = 2.1e5, nu = 0.3;
	fea_model.SetMatrial(em, nu);

	system("pause");
	return 0;
}