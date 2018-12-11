#include <iostream>
// #include <cstdlib>
// #include <vector>
// #include <algorithm>
// #include <mpi.h>
#include <complex>
#include "../common.h"

class SMatrix
{
	// TODO
	int descriptor[9]; // ??
	int myRows, myCols;
	complex_d *data;

	static int rank, nProcs, blockSize, root;
	int icon;
	int nRows, nCols, nProcRows, nProcCols,
		myProcRow, myProcCol,
		myProcRowsOffset, myProcColsOffset,
		myProcRows, myProcCols, myProcSize;
	bool myProcOk;

  public:
	SMatrix();
	SMatrix(int nRows, int nCols, int nProcRows=0, int nProcCols=0);
	SMatrix(int, const char *); // context -- ???
	SMatrix(const SMatrix &);

	~SMatrix();

	SMatrix& operator=(const SMatrix &);
	bool operator ==(const SMatrix&);
	SMatrix& operator +(const SMatrix&);
	SMatrix& operator *(const SMatrix&);
	SMatrix& operator *(complex_d);
	void operator +=(const SMatrix&);
	void operator *=(const SMatrix&);

	// экспонента от матрицы
	static SMatrix& exp(const SMatrix&);
	// сопряжение матрицы
	static SMatrix& conjugate(const SMatrix&);

	// собственные значения, собственные векторы
	void eigen(double *eigenvalues, SMatrix& eigenvectors) const;

	// мб все это внутри?
	bool init(int *rank_p, int *nProcs_p, int blockSize);
	void exit();
	void barrier();

	void read(int context, const char *filename);
	void write(const char *filename);
	void allocate(int context, int N);
	void print_diag() const;
	void print_full() const;
};

void trivial_evolution(int system_size, SMatrix&roStart, SMatrix &H, double dT, int n_steps);
