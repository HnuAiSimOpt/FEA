#ifndef ReadInp_H
#define ReadInp_H

#include <fstream>
#include <iostream>
#include <sstream>
#include"IntMatrix.h"
#include"IntArray.h"
#include"DoubleMatrix.h"


// read total number of elements and nodes  data: 2023/5/10
void ReadNumOfEleNode(string path, int* ne_, int* nd_);

// read node coordinates and node connections  data: 2023/5/10
void ReadInpFile(string path, DoubleMatrix& coordinates, IntMatrix& connections);

// read boundary and load  (data: 2023 / 5 / 11)
void ReadBoundaryLoad(string path, IntArray& constrain, IntArray& load);


#endif
