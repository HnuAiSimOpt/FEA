//
#include"WriteOutputFile.h"


// export displacement  (data: 2023 / 5 / 29)
void ExportDis(string path_, DoubleArray& dis, IntArray& node_map)
{
	string path = path_ + "dis.txt";
	std::ofstream fout(path);
	if (!fout)
	{
		std::cerr << "fail to open file" << std::endl;
	}
	int n_node = node_map.getSize();
	int id = 1;
	for (int i = 1; i <= n_node; i++)
	{
		if (node_map(i))
		{
			fout << dis(3 * id - 2) << "    " << dis(3 * id - 1) << "    " << dis(3 * id) << std::endl;
			id++;
		}
		else
		{
			fout << 0.0 << "    " << 0.0 << "    " << 0.0 << std::endl;
		}
	}
}

// export coordinates  (data: 2023 / 5 / 29)
void ExportDis(string path_, DoubleMatrix& coord)
{
	string path = path_ + "coord.txt";
	std::ofstream fout(path);
	if (!fout)
	{
		std::cerr << "fail to open file" << std::endl;
	}
	int n_node = coord.GetRow();
	int id = 1;
	for (int i = 1; i <= n_node; i++)
	{
		fout << coord(i, 1) << "    " << coord(i, 2) << "    " << coord(i, 3) << std::endl;
	}
}


// export connect  (data: 2023 / 5 / 29)
void ExportConnect(string path_, IntMatrix& connect)
{
	string path = path_ + "connect.txt";
	std::ofstream fout(path);
	if (!fout)
	{
		std::cerr << "fail to open file" << std::endl;
	}
	int n_node = connect.GetRow();
	int id = 1;
	for (int i = 1; i <= n_node; i++)
	{
		fout << connect(i, 1) << "    " << connect(i, 2) << "    " << connect(i, 3) << "    "
			<< connect(i, 4) << "    " << connect(i, 5) << "    " << connect(i, 6) << "    "
			<< connect(i, 7) << "    " << connect(i, 8) << std::endl;
	}
}

// write VYK file (data: 2023 / 5 / 29)
void ExportData2VtkFile(string path_, double factor, DoubleMatrix& coord, IntMatrix& connect, 
	DoubleArray& dis, IntArray& node_map)
{
	string path = path_ + "show.vtk";
	string DATASET[6] = { "STRUCTURED_POINTS", "STRUCTURED_GRID", "UNSTRUCTURED_GRID", "POLYDATA", "RECTILINEAR_GRID", "FIELD" };
	std::ofstream fout;
	fout.open(path, std::ios::out);
	if (!fout) 
	{
		std::cerr << "Error: Cannot open " << path << std::endl;
		exit(EXIT_FAILURE);
	}
	fout << "# vtk DataFile Version 3.0\n";                  // Version Statement
	fout << "The density field of the optimized results\n";  // title
	fout << "ASCII\n";                                       // file format statement
	fout << "DATASET " + DATASET[2] + "\n\n";                // data format: unstructured grid

	// write inform of node coordinates
	int num_node = coord.GetRow();
	fout << "POINTS\t" << num_node << "\tdouble\n";

	int id = 1;
	for (int i = 1; i < num_node + 1; i++) 
	{
		if (node_map(i))
		{
			fout << coord(i, 1) + factor * dis(3 * id - 2) << "\t\t"
				<< coord(i, 2) + factor * dis(3 * id - 1) << "\t\t"
				<< coord(i, 3) + factor * dis(3 * id) << "\n";
			id++;
		}
		else
		{
			fout << coord(i, 1) << "\t\t" << coord(i, 2) << "\t\t" << coord(i, 3) << "\n";
		}
	}

	// write inform of connect relationship
	int num_ele = connect.GetRow();
	fout << "CELLS\t" << num_ele << "\t" << num_ele * (8 + 1) << "\n";
	for (int i = 1; i < num_ele + 1; i++)
	{
		fout << 8 << "\t" << connect(i, 1) - 1 << "\t" << connect(i, 2) - 1 << "\t" << connect(i, 4) - 1 << "\t"
			<< connect(i, 3) - 1 << "\t" << connect(i, 5) - 1 << "\t" << connect(i, 6) - 1 << "\t"
			<< connect(i, 8) - 1 << "\t" << connect(i, 7) - 1 << "\n";
	}
	// write inform of element type
	fout << "CELL_TYPES\t\t" << num_ele << "\n";
	for (int i = 1; i < num_ele + 1; i++)
	{
		fout << 11 << "\n";
	}
	// write infom of node displacement
	fout << "\nPOINT_DATA\t" << num_node << "\nSCALARS ux1 double 1\n" << "LOOKUP_TABLE  table1\n";
	id = 1;
	for (int i = 1; i < num_node + 1; i++)
	{
		if (node_map(i))
		{
			fout << dis(3 * id - 2) << "\n";
			id++;
		}
		else
		{
			fout << 0 << "\n";
		}
	}
	fout << "\nSCALARS uy1 double 1\n" << "LOOKUP_TABLE  table2\n";
	id = 1;
	for (int i = 1; i < num_node + 1; i++)
	{
		if (node_map(i))
		{
			fout << dis(3 * id - 1) << "\n";
			id++;
		}
		else
		{
			fout << 0 << "\n";
		}
	}
	fout << "\nSCALARS uz1 double 1\n" << "LOOKUP_TABLE  table3\n";
	id = 1;
	for (int i = 1; i < num_node + 1; i++)
	{
		if (node_map(i))
		{
			fout << dis(3 * id) << "\n";
			id++;
		}
		else
		{
			fout << 0 << "\n";
		}
	}
}