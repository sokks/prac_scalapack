#ifndef SMATRIX_H
#define SMATRIX_H

#define USE_DOUBLE 1

#include "mpi.h"
#include <ostream>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>

using namespace std;

#ifdef USE_DOUBLE
	typedef double SReal;
	#define pRgesvd pdgesvd_
	#define pCgesvd pzgesvd_
#else
	typedef float SReal;
	#define pRgesvd psgesvd_
	#define pCgesvd pcgesvd_
#endif

#ifdef USE_COMPLEX
	#define pNgesvd pCgesvd
	typedef complex<SReal> SType;
#else
	#define pNgesvd pRgesvd
	typedef SReal SType;
#endif

class SMatrix {
public:
	SMatrix(int nRows, int nCols, int nProcRows=0, int nProcCols=0);
	SMatrix(const SMatrix& matrix);
	SMatrix(const SMatrix& matrix, int nRows, int nCols);
	~SMatrix();

	static bool init(int *rank=NULL, int *nprocs=NULL, int blockSize=64);
	static void exit();
	static void barrier();
	void dumpInfo() const;
	void dump() const;
	void operator>>(ostream& stream) const;

	void setBlock(int procRow, int procCol, const SType *data, int size = -1, int offset = 0);
	void set(const SType *data, int rowStart = 0, int colStart = 0, int rowCount = 0, int colCount = 0);
	void fill(SType value);
	void populate(SType (*func)(int, int));
	void setZero();
	void setIdentity();

	SMatrix& operator =(const SMatrix& matrix);
	bool operator ==(const SMatrix& matrix);
	SMatrix& operator =(SType value);
	SMatrix& operator *=(SType value);
	SMatrix& operator /=(SType value);
	SMatrix& operator +=(SType value);
	SMatrix& operator -=(SType value);

	void readf(MPI_File thefile);
	void wrbuf(SType* buf);

	void setij(int i, int j, SType val);
	void setij(int i, int j);

	void fill_RWA(std::vector<int> H_idxs, std::map<int, int> H_sizes);

private:
	void create(int icon = -2);
	void destroy();
	
	static SType fillerIdentity(int i, int j);
	int* getDesc();
	int* getDesc(int rows, int cols, int procRows);

	static int rank, nProcs, blockSize, root;
	int icon;
	int nRows, nCols, nProcRows, nProcCols,
		myProcRow, myProcCol,
		myProcRowsOffset, myProcColsOffset,
		myProcRows, myProcCols, myProcSize;
	bool myProcOk;
	SType *data;

// for task2
public:
	int N = 0;
	double w_c = 1.0;
	double w_a = 10.0;
	double a = 1.0;
	double b = 100.0;

private:
	std::vector<int> get_H_p_vectors(int p);
	double get_interaction(int vec1, int vec2);
	double **gen_H_p(int p, int sz);
	double get_H_p_i_j(int p, int i, int j, std::map<int, double **>& H_generated, std::map<int, int> H_sizes);

};

ostream& operator<<(ostream& out, const SMatrix& matrix);

#endif
