#pragma once
#include"vtk.h"
#include <iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<iomanip>
using namespace std;

vtk::vtk()
{
	dis_display = false;
	ratio = 0;
}
vtk::vtk(bool dis_display, float ratio)
{
	this->dis_display = dis_display;
	this->ratio = ratio;
}
void vtk::vtkoutput()
{
	string DATASET[6] = { "STRUCTURED_POINTS" ,
							"STRUCTURED_GRID",
							"UNSTRUCTURED_GRID",
							"POLYDATA",
							"RECTILINEAR_GRID",
							"FIELD"
	};
	ifstream coord, dis, connect, stress;
	ofstream fout;
	fout.open("show.vtk", ios::out);
	if (!fout)
	{
		cout << "����ļ���ʧ��" << endl;

	}
	fout << "# vtk DataFile Version 3.0\n"; //�汾����
	fout << "The density field of the optimized results\n";  //����
	fout << "ASCII\n";   //�ļ���ʽ����
	fout << "DATASET " + DATASET[2] + "\n\n";    //���ݸ�ʽ���ǽṹ������
	coord.open("coord.txt", ios::in);
	dis.open("dis.txt", ios::in);
	connect.open("connect.txt", ios::in);
	//stress.open("stress.txt", ios::in);
	if (!coord)
	{
		cout << "coord.txt�ļ���ʧ�ܣ�" << endl;

	}
	if (!dis)
	{
		cout << "dis.txt�ļ���ʧ�ܣ�" << endl;

	}
	if (!connect)
	{
		cout << "connect.txt�ļ���ʧ�ܣ�" << endl;

	}
	vector<string> out_line;
	string line;
	string data;
	string num_coor;
	string num_dis;
	string num_ele;
	string num_stress;
	// ��ȡcoor�ļ�   �ڵ�������Ϣ
	if (this->dis_display)
	{
		string distance;
		vector<string> out_dis;
		double x_new, y_new, z_new;
		while (!coord.eof())
		{
			getline(dis, distance);
			getline(coord, line);
			if (line != "")
			{
				stringstream s;
				stringstream s1;
				s << line;
				s1 << distance;
				while (s >> data)
				{
					out_line.push_back(data);
				}
				while (s1 >> data)
				{
					out_dis.push_back(data);
				}
				if (out_line[0] == "POINTS")
				{
					num_coor = out_line[1];
					fout << "POINTS\t" + num_coor + '\t' + "double\n";
				}
				else
				{
					x_new = atof(out_line[0].c_str()) + atof(out_dis[0].c_str()) * this->ratio;
					y_new = atof(out_line[1].c_str()) + atof(out_dis[1].c_str()) * this->ratio;
					z_new = atof(out_line[2].c_str()) + atof(out_dis[2].c_str()) * this->ratio;
					fout << std::setw(20) << std::setiosflags(ios::left) << setprecision(4) << x_new << std::setw(20)
						<< y_new << std::setw(20) << z_new << "\n";
				}
				out_line.clear();
				out_dis.clear();
			}
		}
		dis.clear();
		dis.close();
		dis.open("dis.txt", ios::in);
	}
	else
	{
		while (!coord.eof())
		{
			getline(coord, line);
			if (line != "")
			{
				stringstream s;
				s << line;
				while (s >> data)
				{
					out_line.push_back(data);
				}
				if (out_line[0] == "POINTS")
				{
					num_coor = out_line[1];
					fout << "POINTS\t" + num_coor + '\t' + "double\n";
				}
				else
				{
					fout << std::setw(20) << std::setiosflags(ios::left) << setprecision(4) << out_line[0] << std::setw(20)
						<< out_line[1] << std::setw(20) << out_line[2] << "\n";
				}
				out_line.clear();
			}
		}
	}

	// ��ȡconnect�ļ�    ��Ԫ�ڵ���Ϣ
	int Elenode = 8;	//��Ԫ�ڵ����
	int Elenum;//��Ԫ����
	fout << "\n";
	while (!connect.eof())
	{
		getline(connect, line);
		if (line != "")
		{
			stringstream s;
			s << line;
			while (s >> data)
			{
				out_line.push_back(data);
			}
			if (out_line[0] == "CELLS")
			{
				num_ele = out_line[1];
				Elenum = atof(num_ele.c_str());
				fout << "CELLS\t" + num_ele + '\t' << Elenum * (1 + Elenode) << "\n";
			}
			else
			{
				fout << std::setw(10) << Elenode;
				fout << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[0].c_str()) - 1 << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[1].c_str()) - 1;
				fout << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[3].c_str()) - 1 << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[2].c_str()) - 1;
				fout << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[4].c_str()) - 1 << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[5].c_str()) - 1;
				fout << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[7].c_str()) - 1 << std::setw(10) << std::setiosflags(ios::left) << atof(out_line[6].c_str()) - 1;
				fout << "\n";
			}
			out_line.clear();
		}
	}
	//�����Ԫ����
	fout << "\n";
	fout << "CELL_TYPES\t" << num_ele << "\n";
	string type;
	switch (Elenode)
	{
	case 8:	//������
	{
		type = "11";
		break;
	}
	case 4:	//������
	{
		type = "10";
		break;
	}

	default:
		break;
	}
	for (int i = 0; i < Elenum; i++)
	{
		fout << type << "\n";
	}

	// ��ȡdis�ļ�    �ڵ�λ����Ϣ
	fout << "\n";
	fout << "POINT_DATA\t" << num_coor << "\n";
	fout << "SCALARS ux1 double 1\n";
	fout << "LOOKUP_TABLE  table1\n";
	while (!dis.eof())
	{
		getline(dis, line);
		if (line != "")
		{
			stringstream s;
			s << line;
			while (s >> data)
			{
				out_line.push_back(data);
			}
			if (out_line[0] == "POINTS")
			{
				num_dis = out_line[1];
				if (num_dis != num_coor)
				{
					cout << "coor�ļ���dis�ļ���POINTSֵ��ͬ" << endl;

				}
			}
			else
			{
				fout << std::setiosflags(ios::left) << out_line[0] << "\n";
			}
			out_line.clear();
		}
	}
	dis.clear();
	dis.close();
	dis.open("dis.txt", ios::in);
	if (!dis.is_open())
	{
		cout << "dis.txt�ڶ������´�ʧ��" << endl;
	}
	fout << "\n";
	fout << "SCALARS uy1 double 1\n";
	fout << "LOOKUP_TABLE  table2\n";
	while (!dis.eof())
	{
		getline(dis, line);
		if (line != "")
		{
			stringstream s;
			s << line;
			while (s >> data)
			{
				out_line.push_back(data);
			}
			if (out_line[0] == "POINTS")
			{
				out_line.clear();
				continue;
			}
			else
			{
				fout << std::setiosflags(ios::left) << out_line[1] << "\n";
			}
			out_line.clear();
		}
	}
	dis.clear();
	dis.close();
	dis.open("dis.txt", ios::in);
	if (!dis.is_open())
	{
		cout << "dis.txt���������´�ʧ��" << endl;
	}
	fout << "\n";
	fout << "SCALARS uz1 double 1\n";
	fout << "LOOKUP_TABLE  table3\n";
	while (!dis.eof())
	{
		getline(dis, line);
		if (line != "")
		{
			stringstream s;
			s << line;
			while (s >> data)
			{
				out_line.push_back(data);
			}
			if (out_line[0] == "POINTS")
			{
				out_line.clear();
				continue;
			}
			else
			{
				fout << std::setiosflags(ios::left) << out_line[2] << "\n";
			}
			out_line.clear();
		}
	}
	//��ȡstress�ļ�    ��ԪӦ����Ϣ
	fout << "CELL_DATA\t" << num_ele << "\n";
	fout << "SCALARS stress double 1\n";
	fout << "LOOKUP_TABLE  table 1\n";
	while (!stress.eof())
	{
		getline(stress, line);
		if (line != "")
		{
			stringstream s;
			s << line;
			while (s >> data)
			{
				out_line.push_back(data);
			}
			if (out_line[0] == "CELLS")
			{
				num_stress = out_line[1];
				if (num_stress != num_ele)
				{
					cout << "stress�ļ���connect�ļ���CELLֵ��ͬ" << endl;

				}
			}
			else
			{
				fout << std::setiosflags(ios::left) << out_line[0] << "\n";
			}
			out_line.clear();
		}
	}
	fout.close();
	coord.close();
	connect.close();
	dis.close();

}