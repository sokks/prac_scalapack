#include <iostream>
// #include <cstdlib>
// #include <vector>
// #include <algorithm>
// #include <mpi.h>
#include <complex>
#include "../common.h"

class SMatrix
{

	int descriptor[9]; // ??
	int myRows, myCols;
	complex_d *myData;

  public:
	SMatrix();
	SMatrix(int, const char *); // context -- ???
	SMatrix(const SMatrix &);

	~SMatrix();

	SMatrix &operator=(const SMatrix &A);

	void read(int context, const char *filename);
	void write(const char *filename);
	void allocate(int context, int N);
	void print_diag();
	void print_full();
};
