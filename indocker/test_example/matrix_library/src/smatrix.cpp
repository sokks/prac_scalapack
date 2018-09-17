#include "smatrix.h"
#include "scalapack.h"

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

int SMatrix::rank, SMatrix::nProcs, SMatrix::blockSize, SMatrix::root = 0;

SMatrix::SMatrix(int nRows, int nCols, int nProcRows, int nProcCols)
: nRows(nRows), nCols(nCols)
{
	if ((nProcRows > nProcs) || (nProcCols > nProcs)) {
		cout << "Too many processes requested." << endl;
		::exit(1);
	}
	if ((nProcRows <= 0) && (nProcCols <= 0)) {
		int divA;
		for (divA = (int) sqrt(nProcs); (nProcs % divA != 0); divA--);
		int divB = nProcs / divA;
		cout << divA << endl;
		if (nRows > nCols) {
			nProcRows = divB;
			nProcCols = divA;
		} else {
			nProcRows = divA;
			nProcCols = divB;
		}
	} else if (nProcRows <= 0) {
		nProcRows = nProcs / nProcCols;
	} else if (nProcCols <= 0) {
		nProcCols = nProcs / nProcRows;
	} else if (nProcs < nProcRows * nProcCols) {
		cout << "Too many processes requested." << endl;
		::exit(1);
	}

	this->nProcRows = nProcRows;
	this->nProcCols = nProcCols;

	create();
}

SMatrix::SMatrix(const SMatrix &matrix)
: nRows(matrix.nRows), nCols(matrix.nCols), nProcRows(matrix.nProcRows), nProcCols(matrix.nProcCols)
{
	create(matrix.icon);
	memcpy(data, matrix.data, myProcSize * sizeof(SType));
}

SMatrix::SMatrix(const SMatrix &matrix, int nRows, int nCols)
: nRows(nRows), nCols(nCols), nProcRows(matrix.nProcRows), nProcCols(matrix.nProcCols)
{
	create(matrix.icon);
}

SMatrix::~SMatrix() {
	destroy();
}

void SMatrix::create(int icon) {
	if (icon < -1) {
		Cblacs_get(-1, 0, &icon);
		Cblacs_gridinit(&icon, (char *) "Column", nProcRows, nProcCols);
	}
	this->icon = icon;
	int  nProcRows_, nProcCols_;
	Cblacs_gridinfo(icon, &nProcRows_, &nProcCols_, &myProcRow, &myProcCol);
	myProcOk = (myProcRow > -1) && (myProcCol > -1) && (myProcRow < nProcRows_) & (myProcCol < nProcCols_);
	if (myProcOk) {
		myProcRows = numroc_(&nRows, &blockSize, &myProcRow, &root, &nProcRows);
		myProcCols = numroc_(&nCols, &blockSize, &myProcCol, &root, &nProcCols);
		myProcRowsOffset = npreroc_(&nRows, &blockSize, &myProcRow, &root, &nProcRows);
		myProcColsOffset = npreroc_(&nCols, &blockSize, &myProcCol, &root, &nProcRows);
	} else {
		myProcRows = myProcCols = myProcRowsOffset = myProcColsOffset = 0;
	}
	myProcSize = myProcRows * myProcCols;
	if (myProcSize == 0) {
		data = NULL;
		return;
	}
	data = new SType[myProcSize];
}
void SMatrix::destroy() {
	delete[] data;
}

void SMatrix::setBlock(int procRow, int procCol, const SType *data, int size, int offset) {
	if ((procRow != myProcRow) || (procCol != myProcCol)) {
		return;
	}
	if (size == -1) {
		size = myProcSize - offset;
	}
	memcpy(this->data + offset, &data, size * sizeof(SType));
}

void SMatrix::set(const SType *data, int rowStart, int colStart, int rowCount, int colCount) {
	if (rowCount == 0) {
		rowCount = nRows - rowStart;
	}
	if (colCount == 0) {
		colCount = nCols - colStart;
	}
	int	rowEnd = rowStart + rowCount,
		colEnd = colStart + colCount;
	if (	(rowEnd <= myProcRowsOffset) || 
		(colEnd <= myProcColsOffset) ||
		(rowStart >= myProcRowsOffset + myProcRows) ||
		(colStart >= myProcColsOffset + myProcCols)) {
		return;
	}

	for (int row = max(myProcRowsOffset, rowStart), myRow = row - myProcRowsOffset; (myRow < myProcRows) && (row < rowEnd); myRow++, row++) {
		// memcpy won't work here, because of different-major order
		for (int col = max(myProcColsOffset, colStart), myCol = col - myProcColsOffset; (myCol < myProcCols) && (col < colEnd); myCol++, col++) {
			this->data[myRow + myProcRows * myCol] = data[colCount * (row - rowStart) + (col - colStart)];
		}
	}
}

void SMatrix::fill(SType value) {
	barrier();
	for (int i = 0; i < myProcSize; i++) {
		data[i] = value;
	}
	barrier();
}

void SMatrix::populate(SType (*func)(int, int)) {
	barrier();
	for (int i = 0; i < myProcRows; i++) {
		for (int j = 0; j < myProcCols; j++) {
			data[i + myProcRows * j] = func(i + myProcRowsOffset, j + myProcColsOffset);
		}
	}
	barrier();
}

void SMatrix::setZero() {
	fill(0);
}

void SMatrix::setIdentity() {
	populate(&fillerIdentity);
}

SMatrix& SMatrix::operator =(const SMatrix& matrix) {
	destroy();
	nRows = matrix.nRows;
	nCols = matrix.nCols;
	nProcRows = matrix.nProcRows;
	nProcCols = matrix.nProcCols;
	create(matrix.icon);
	memcpy(data, matrix.data, myProcSize * sizeof(SType));
}
bool SMatrix::operator ==(const SMatrix& matrix) {
	unsigned char result, equals = (memcmp(data, matrix.data, myProcSize * sizeof(SType)) == 0);
	MPI_Allreduce(&equals, &result, 1, MPI_UNSIGNED_CHAR, MPI_LAND, MPI_COMM_WORLD);
	return (result != 0);
}

SMatrix& SMatrix::operator =(SType value) {
	fill(value);
	return *this;
}

SMatrix& SMatrix::operator *=(SType value) {
	barrier();
	for (int i = 0; i < myProcSize; i++) {
		data[i] *= value;
	}
	barrier();
	return *this;
}

SMatrix& SMatrix::operator /=(SType value) {
	barrier();
	for (int i = 0; i < myProcSize; i++) {
		data[i] /= value;
	}
	barrier();
	return *this;
}

SMatrix& SMatrix::operator +=(SType value) {
	barrier();
	for (int i = 0; i < myProcSize; i++) {
		data[i] += value;
	}
	barrier();
	return *this;
}

SMatrix& SMatrix::operator -=(SType value) {
	barrier();
	for (int i = 0; i < myProcSize; i++) {
		data[i] -= value;
	}
	barrier();
	return *this;
}

SType SMatrix::fillerIdentity(int i, int j) {
	return (i == j) ? 1 : 0;
}

int* SMatrix::getDesc() {
	return getDesc(nRows, nCols, myProcRows);
}

int* SMatrix::getDesc(int rows, int cols, int procRows) {
	int *desc = new int[9];
	int itemp = max(1, procRows), izero = 0, info;
	descinit_(desc, &rows, &cols, &blockSize, &blockSize, &izero, &izero, &icon, &itemp, &info);
	return desc;
}

SReal* SMatrix::calculateSVD(int *SSize, SMatrix **U, SMatrix **VT) {
	int nMin = min(nRows, nCols);
	SMatrix	*u = new SMatrix(*this, nRows, nMin),
		*vt = new SMatrix(*this, nMin, nCols);
	SReal *S = new SReal[nMin];

        int	*desc = getDesc(),
		*descU = u->getDesc(),
		*descVT = vt->getDesc();

	SType *dataCopy = new SType[myProcSize];
	memcpy(dataCopy, data, myProcSize * sizeof(SType));

	SType *work = new SType[1];

	int ione = 1, info, lwork = -1;

#ifdef USE_COMPLEX
	SReal *rwork = new SReal[1 + 4 * max(nRows, nCols)];
#endif
	pNgesvd((char *) "V", (char *) "V", &nRows, &nCols, dataCopy, &ione, &ione, desc,
		S, u->data, &ione, &ione, descU, vt->data, &ione, &ione, descVT,
		work, &lwork,
#ifdef USE_COMPLEX
		rwork,
#endif
		&info);

#ifdef USE_COMPLEX
	lwork = (int) work[0].real();
#else
	lwork = (int) work[0];
#endif

	delete[] work;
	work = new SType[lwork];

	pNgesvd((char *) "V", (char *) "V", &nRows, &nCols, dataCopy, &ione, &ione, desc,
		S, u->data, &ione, &ione, descU, vt->data, &ione, &ione, descVT,
		work, &lwork,
#ifdef USE_COMPLEX
		rwork,
#endif
		&info);

	delete[] work, dataCopy;
#ifdef USE_COMPLEX
	delete[] rwork;
#endif
	delete[] desc, descU, descVT;

	if (SSize != NULL) {
		*SSize = nMin;
	}

	if (U != NULL) {
		*U = u;
	} else {
		delete u;
	}

	if (VT != NULL) {
		*VT = vt;
	} else {
		delete vt;
	}

	return S;
}
  
bool SMatrix::init(int *rank_p, int *nProcs_p, int blockSize) {
	SMatrix::blockSize = blockSize;
	Cblacs_pinfo(&rank, &nProcs);
	if (rank_p != 0) {
		*rank_p = rank;
	}
	if (nProcs_p != 0) {
		*nProcs_p = nProcs;
	}
	return (rank == root);
}

void SMatrix::exit() {
//	Cblacs_gridexit(0);
	Cblacs_exit(0);
}

void SMatrix::barrier() {
	MPI_Barrier(MPI_COMM_WORLD);
}


void SMatrix::dumpInfo() const {
	cout << flush;
	barrier();
	if (rank == root) {
		cout << endl;
		cout << "context: " << icon << ", blockSize: " << blockSize << endl;
		cout << "nRows: " << nRows << ", nCols: " << nCols << endl;
		cout << "nProcRows: " << nProcRows << ", nProcCols: " << nProcCols << endl;
		cout << endl;
		cout << "procId\trows\tcols\trowFrom\tcolFrom\tsize" << endl;
		cout << flush;
	}
	for (int i=0; i<nProcs; i++) {
		if (rank == i) {
			cout	<< rank << "\t" << myProcRows << "\t" << myProcCols << "\t"
				<< myProcRowsOffset << "\t" << myProcColsOffset << "\t" << myProcSize
				<< endl << flush;
		}
		barrier();
	}
	if (rank == root) {
		cout << endl << flush;
	}
}
void SMatrix::dump() const {
	cout << flush;
	barrier();
	if (rank == root) {
		cout << "nProcRows=" << nProcRows << ", nProcCols=" << nProcCols << endl << flush;
	}
	for (int i=0; i<nProcRows; i++) {
		for (int j=0; j<nProcCols; j++) {
			if ((i == myProcRow) && (j == myProcCol)) {
				cout << "block(" << i << ", " << j << "):";
				for (int k = 0; k < myProcSize; k++) {
					cout
#ifdef USE_COMPLEX
						<< setprecision(2) << setw(12)
#else
						<< setprecision(4) << setw(9)
#endif
						<< data[k];
				}
				cout << endl << flush;
			}
			barrier();
		}
	}
}
void SMatrix::operator>>(ostream& out) const {
	cout << flush;
	barrier();
	for (int i=0; i<nRows; i++) {
		for (int j=0; j<nCols; j++) {
			// Can use indxg2p_ and indxg2l_ here, but do not want to.
			if (	(i >= myProcRowsOffset) &&
				(i < myProcRowsOffset + myProcRows) &&
				(j >= myProcColsOffset) &&
				(j < myProcColsOffset + myProcCols)) {
				cout
#ifdef USE_COMPLEX
					<< setprecision(2) << setw(12)
#else
					<< setprecision(4) << setw(9)
#endif

					<< data[(i - myProcRowsOffset) + myProcRows * (j - myProcColsOffset)];
				if (j == (nCols - 1)) {
					cout << endl;
				}
				cout << flush;
			}
			barrier();
		}
	}
}
ostream& operator<<(ostream& out, const SMatrix& matrix) {
	matrix >> out;
	return out;
}
