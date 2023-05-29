#ifndef WRITEOUTPUTFILE_H
#define WRITEOUTPUTFILE_H

#include<string>
#include<fstream>
#include<vector>
#include<iomanip>
#include <iostream>
#include"IntArray.h"
#include"IntMatrix.h"
#include"DoubleArray.h"
#include"DoubleMatrix.h"

// export displacement  (data: 2023 / 5 / 29)
void ExportDis(string path, DoubleArray& dis, IntArray& reset_idx);

// export coordinates  (data: 2023 / 5 / 29)
void ExportDis(string path, DoubleMatrix& coord);

// export connect  (data: 2023 / 5 / 29)
void ExportConnect(string path, IntMatrix& connect);

// write VYK file (data: 2023 / 5 / 29)
void ExportData2VtkFile(string path_, double factor, DoubleMatrix& coord, IntMatrix& connect,
	DoubleArray& dis, IntArray& reset_idx);

#endif
