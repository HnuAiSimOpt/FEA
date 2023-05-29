#include"ReadInp.h"

// read total number of elements and nodes  data: 2023/5/10
void ReadNumOfEleNode(string path, int* ne_, int* nd_)
{
	std::ifstream infile(path.c_str(), std::ios::in);
	if (!infile)
	{
		std::cerr << "Error: Cannot open " << path << std::endl;
		exit(EXIT_FAILURE);
	}
	// read total element numbers and total numbers
	string line;
	while (getline(infile, line))
	{
		if (line.find("*Nset") != string::npos)
		{
			getline(infile, line);
			std::istringstream iss(line);
			string t1, nd, t2;
			iss >> t1 >> nd >> t2;
			*nd_ = atoi(nd.c_str());
		}
		if (line.find("*Elset") != string::npos)
		{
			getline(infile, line);
			std::istringstream iss(line);
			string t1, ne, t2;
			iss >> t1 >> ne >> t2;
			*ne_ = atoi(ne.c_str());
			break;
		}
	}
	infile.close();
}


// read node coordinates and node connections data: 2023/5/10
void ReadInpFile(string path, DoubleMatrix& coordinates, IntMatrix& connections)
{
	// read inp file
	std::ifstream infile(path.c_str(), std::ios::in);  
	string line;
	// record data
	while (getline(infile, line))
	{
		if (line.find("*Node") != string::npos)
		{
			int id_node = 1;
			while (getline(infile, line))
			{
				if (line.find("*") != string::npos)
					break;
				else
				{
					string id, x, y, z;
					double x_, y_, z_;
					std::istringstream iss(line);
					iss >> id >> x >> y >> z;
					x.erase(x.end() - 1);  // delete the last symbol in the string
					x_ = stod(x);          // convert to double type
					y.erase(y.end() - 1);
					y_ = stod(y);
					z_ = stod(z);          // no symbol in the string
					coordinates.SetValues(id_node, 1, x_);
					coordinates.SetValues(id_node, 2, y_);
					coordinates.SetValues(id_node, 3, z_);
					id_node = id_node + 1;
				}
			}
		}
		if (line.find("*Element") != string::npos)
		{
			int id_ele = 1;
			while (getline(infile, line))
			{
				if (line.find("*") != string::npos)
					break;
				else
				{
					string id, node1, node2, node3, node4, node5, node6, node7, node8;
					int node1_, node2_, node3_, node4_, node5_, node6_, node7_, node8_;
					std::istringstream iss(line);
					iss >> id >> node1 >> node2 >> node3 >> node4 >> node5 >> node6 >> node7 >> node8;
					node1.erase(node1.end() - 1);  // delete the last symbol in the string
					node1_ = atoi(node1.c_str());  // convert to int type
					node2.erase(node2.end() - 1);
					node2_ = atoi(node2.c_str());
					node3.erase(node3.end() - 1);
					node3_ = atoi(node3.c_str());
					node4.erase(node4.end() - 1);
					node4_ = atoi(node4.c_str());
					node5.erase(node5.end() - 1);
					node5_ = atoi(node5.c_str());
					node6.erase(node6.end() - 1);
					node6_ = atoi(node6.c_str());
					node7.erase(node7.end() - 1);
					node7_ = atoi(node7.c_str());
					node8_ = atoi(node8.c_str());
					connections.SetValues(id_ele, 1, node1_);
					connections.SetValues(id_ele, 2, node2_);
					connections.SetValues(id_ele, 3, node3_);
					connections.SetValues(id_ele, 4, node4_);
					connections.SetValues(id_ele, 5, node5_);
					connections.SetValues(id_ele, 6, node6_);
					connections.SetValues(id_ele, 7, node7_);
					connections.SetValues(id_ele, 8, node8_);
					id_ele = id_ele + 1;
				}
			}
		}
		if (line.find("*End Instance") != string::npos)
		{
			break;
		}
	}
	infile.close();
}

// read boundary and load  (data: 2023 / 5 / 11)
void ReadBoundaryLoad(string path, IntArray& constrain, IntArray& load)
{
	// read inp file
	std::ifstream infile(path.c_str(), std::ios::in);
	string line;
	// record data
	while (getline(infile, line))
	{
		string node_id;
		int node_id_;
		if (line.find("Set-constrain") != string::npos)
		{
			while (getline(infile, line))
			{
				if (line.find("*") != string::npos)
					break;
				else
				{
					std::istringstream iss(line);
					while (iss >> node_id)
					{
						if (node_id.find(",") != string::npos)
						{
							node_id.erase(node_id.end() - 1);   // delete the last symbol in the string
						}
						node_id_ = atoi(node_id.c_str());       // convert to int type
						constrain.push_back(node_id_);
					}
				}
			}
		}
		if (line.find("Set-load") != string::npos)
		{
			while (getline(infile, line))
			{
				if (line.find("*") != string::npos)
					break;
				else
				{
					std::istringstream iss(line);
					while (iss >> node_id)
					{
						if (node_id.find(",") != string::npos)
						{
							node_id.erase(node_id.end() - 1);   // delete the last symbol in the string
						}
						node_id_ = atoi(node_id.c_str());       // convert to int type
						load.push_back(node_id_);
					}
				}
			}
		}
		if (line.find("*End Assembly") != string::npos)
		{
			break;
		}
	}
	infile.close();
	constrain.sortmin2max();
	load.sortmin2max();
}