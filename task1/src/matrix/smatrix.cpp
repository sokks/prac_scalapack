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
		//cout << divA << endl;
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
		myProcColsOffset = npreroc_(&nCols, &blockSize, &myProcCol, &root, &nProcCols);
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
	if (data != NULL) {
		delete[] data;
		data = NULL;
	}
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
	myProcRow = matrix.myProcRow;
	myProcCol = matrix.myProcCol;
	myProcColsOffset = matrix.myProcColsOffset;
	myProcRowsOffset = matrix.myProcRowsOffset;
	myProcRows = matrix.myProcRows;
	myProcCols = matrix.myProcCols;
	myProcSize = matrix.myProcSize;
	create(matrix.icon);
	memcpy(data, matrix.data, myProcSize * sizeof(SType));
	return *this;
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
	if (info != 0) {
		std::cout << "ERROR: descinit_ RETURN CODE: " << info << std::endl;
	} 
	// else {
	// 	std::cout << "OK " << info << std::endl;
	// }
	return desc;
}

void SMatrix::readf(MPI_File thefile) {
	MPI_Status status;
	SType data[myProcSize];
	MPI_Datatype newtype1, newtype;
	MPI_Type_vector(myProcRows, 1, nCols, MPI_DOUBLE_COMPLEX, &newtype1);
	MPI_Type_commit(&newtype1);
	MPI_Type_hvector(myProcCols, 1, sizeof(SType), newtype1, &newtype);
	MPI_Type_commit(&newtype);
	MPI_File_set_view(thefile, 2*sizeof(int) + (nCols*myProcRowsOffset+myProcColsOffset)*sizeof(SType),MPI_DOUBLE_COMPLEX, newtype, "native", MPI_INFO_NULL);
	MPI_File_read(thefile, data, myProcCols*myProcRows, MPI_DOUBLE_COMPLEX, &status);
	barrier();
	for (int i = 0; i < myProcSize; i++) {
		this->data[i] = data[i];
	}
	barrier();
    return;
}

void SMatrix::wrbuf(SType* buf) {
	for (int i = 0; i < nRows; i++) {
			if ((i >= myProcRowsOffset) &&
				(i < myProcRowsOffset + myProcRows) &&
				(i >= myProcColsOffset) &&
				(i < myProcColsOffset + myProcCols)){
					buf[i] = data[(i - myProcRowsOffset) + myProcRows * (i - myProcColsOffset)];
					//cout<<"rank "<<rank<<" i,j "<<i<<" data: "<<buf[i]<<'\n';
				}
	}
    return;
}

void SMatrix::setij(int i, int j) {
	if ((i >= myProcRowsOffset) &&
		(i < myProcRowsOffset + myProcRows) &&
		(j >= myProcColsOffset) &&
		(j < myProcColsOffset + myProcCols)){
		data[(i - myProcRowsOffset) + myProcRows * (j - myProcColsOffset)] = 1;
				}
    return;
}

void SMatrix::mul(char *ch, char *ch_B, SMatrix *B, SMatrix **C) {
	SMatrix	*c = new SMatrix(*this, nRows, B->nCols);
    
	int	*desc  = getDesc();
	int *descB = B->getDesc();
	int	*descC = c->getDesc();

	char *transa = ch;
	char *transb = ch_B;
	
	int m = nRows;
	int n = nCols;
	int k = B->nCols;

	double alpha = 1;
	double beta  = 0;
	int    ione  = 1;


	//cout<<endl<<nRows<<(B)->nCols<<nCols;
	//SType *dataCopy = new SType[myProcSize];
	//memcpy(dataCopy, data, myProcSize * sizeof(SType));
	barrier();
#ifdef USE_COMPLEX
/*
	void pzgemm_ (char *transa, char *transb, int *m, int *n, int *k, 
			double *alpha, 	complex_d *a, int *ia, int *ja, int *desca, 
							complex_d *b, int *ib, int *jb, int *descb, double *beta, 
							complex_d *c, int *ic, int *jc, int *descc);
*/
	pzgemm_(transa, transb, &m, &n, &k, 
			&alpha, 
			data,    &ione, &ione, desc, 
			B->data, &ione, &ione, descB, 
			&beta, 
			c->data, &ione, &ione, descC);
#else
	pdgemm_ ((char *) "N", (char *) "N", &nRows, &(*B)->nCols, &nCols, &alpha, data, &ione, &ione,
		desc, (*B)->data, &ione, &ione, descB, &beta, (*C)->data, &ione, &ione, descC);
#endif
	//delete[] dataCopy;
	delete[] desc, descB, descC;
	// cout<<*c<<endl;
	if (C != NULL) {
		*C = c;
	} else {
		delete c;
	}
	return;
}

SReal* SMatrix::calculateEIGENVAL(int *SSize, SMatrix **Z) {
	//int nMin = min(nRows, nCols);
	SMatrix	*z = new SMatrix(*this, nRows, nRows);
	SReal *W = new SReal[nRows];

        int	*desc = getDesc(),
		*descZ = z->getDesc();

	SType *dataCopy = new SType[myProcSize];
	memcpy(dataCopy, data, myProcSize * sizeof(SType));

	SType *work = new SType[1];

	int ione = 1, info, lwork = -1, lrwork = -1, liwork = -1;

	SReal *rwork = new SReal[1];

	int *iwork = new int[1];
	pzheevd_((char *) "V", (char *) "U", &nRows, dataCopy, &ione, &ione, desc,
		W, z->data, &ione, &ione, descZ, work, &lwork, rwork, &lrwork,
		iwork, &liwork, &info);

#ifdef USE_COMPLEX
	lwork = (int) work[0].real();
#else
	lwork = (int) work[0];
#endif

	delete[] work;
	work = new SType[lwork];
	
	lrwork = (int) rwork[0];
	delete[] rwork;
	rwork = new SReal[lrwork];

	liwork = (int) iwork[0];
	delete[] iwork;
	iwork = new int[liwork];

	pzheevd_((char *) "V", (char *) "U", &nRows, dataCopy, &ione, &ione, desc,
		W, z->data, &ione, &ione, descZ, work, &lwork, rwork, &lrwork,
		iwork, &liwork, &info);

	delete[] work, rwork, iwork, dataCopy;
	delete[] desc, descZ;

	if (SSize != NULL) {
		*SSize = nRows;
	}
	
	if (Z != NULL) {
		*Z = z;
	} else {
		delete z;
	}

	return W;
}

void SMatrix::transpconj(bool isRoot) {

	MPI_Status status;
	SType B[myProcSize];

	if (rank/nProcRows != (rank % nProcRows)){	
		#ifdef USE_COMPLEX
		MPI_Sendrecv(data, myProcSize, MPI_DOUBLE_COMPLEX, rank/nProcRows + (rank % nProcRows) * nProcRows, 1, B, myProcSize, 				MPI_DOUBLE_COMPLEX, rank/nProcRows + (rank % nProcRows) * nProcRows, 1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		#else
		MPI_Sendrecv(data, myProcSize, MPI_DOUBLE, rank/nProcRows + (rank % nProcRows) * nProcRows, 1, B, myProcSize, 				MPI_DOUBLE, rank/nProcRows + (rank % nProcRows) * nProcRows, 1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		#endif
	}
	else {
		for (int i=0; i < myProcSize; i++) {
			B[i] = data[i];
		}
	}

	for (int i=0; i < myProcRows; i++) {
		for (int j=0; j < myProcCols; j++) {
			#ifdef USE_COMPLEX
			//data[j + i*myProcCols] = conj(B[j*myProcRows + i]);
			data[j*myProcRows + i] = conj(B[j + i*myProcCols]);
			#else
			//data[j + i*myProcCols] = B[j*myProcRows + i];
			#endif
		}
	}
	return;
}

void SMatrix::DiagToMat(int *SSize, SType *D){
	barrier();
	(*this).setZero();
	for (int i=0; i<nRows; i++) {
		for (int j=0; j<nCols; j++) {
			// Can use indxg2p_ and indxg2l_ here, but do not want to.
			if (	(i >= myProcRowsOffset) &&
				(i < myProcRowsOffset + myProcRows) &&
				(j >= myProcColsOffset) &&
				(j < myProcColsOffset + myProcCols) &&
				(i==j) ) {
					//x.data[(i - myProcRowsOffset) + myProcRows * (j - myProcColsOffset)] = D[i];
					data[(i - myProcRowsOffset) + myProcRows * (j - myProcColsOffset)] = D[i];
			}
		}
	}
	barrier();
	return;
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
