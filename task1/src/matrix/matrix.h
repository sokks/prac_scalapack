#include <iostream>
// #include <cstdlib>
// #include <vector>
// #include <algorithm>
// #include <mpi.h>
#include <complex>
#include "common.h"

class Matrix {

	int         descriptor[9]; // ??
	int myRows, myCols;
    complex_d * myData;
	
public:
	Matrix();
	Matrix(int, const char *); // context -- ???
	Matrix(Matrix &);

	void read(int context, const char *filename);
	void write(const char *filename);
	void allocate(int context, int N);
	void print();
	void print_full();

	void operator= (const Matrix &A);

	~Matrix();
};
