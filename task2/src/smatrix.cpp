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

void SMatrix::setij(int i, int j, SType val) {
	if ((i >= myProcRowsOffset) &&
		(i < myProcRowsOffset + myProcRows) &&
		(j >= myProcColsOffset) &&
		(j < myProcColsOffset + myProcCols)){
		data[(i - myProcRowsOffset) + myProcRows * (j - myProcColsOffset)] = 1;
				}
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





int get_H_idx(int global_i, int global_j, std::vector<int> H_idxs, std::map<int, int> H_sizes, int *H_i, int *H_j) {
    int idx_i = -1;
    int idx_j = -1;
    
    int prev_sum_i = 0;
    for (int i = 0; i < H_idxs.size(); i++) {
        int sz = H_sizes[H_idxs[i]];
        if (global_i < prev_sum_i + sz) {
            idx_i = H_idxs[i];
            break;
        }
        prev_sum_i += sz;
    }

    int prev_sum_j = 0;
    for (int j = 0; j < H_idxs.size(); j++) {
        int sz = H_sizes[H_idxs[j]];
        if (global_j < prev_sum_j + sz) {
            idx_j = H_idxs[j];
            break;
        }
        prev_sum_j += sz;
    }

    if ( (idx_i != -1) && (idx_j != -1) && (idx_i == idx_j) ) {
        *H_i = global_i - prev_sum_i;
        *H_j = global_j - prev_sum_j;
        return idx_i;
    } 

    return -1;
}


int get_number_of_ones(int n) {
    unsigned int count = 0;
    int n_start = n;
    for (; n; n >>= 1) {
        count += n & 1;
    }

    return count;
}

int get_number_of_ones_odd(int n) {
    unsigned int count = 0;
    for (n = n >> 1; n; n >>= 2) {
        count += n & 1;
    }
    return count;
}

// количество единиц на четных местах (если считать с конца с 0) 
// н-р: [00...00101 -> 2] [00...00010 -> 0]
int get_number_of_ones_even(int n) {
    unsigned int count = 0;
    for (; n; n >>= 2) {
        count += n & 1;
    }
    return count;
}

int power_2_n(int n) {
    int res = 1;
    for (int i = 0; i < n; i++) {
        res <<= 1;
    }
    return res;
}

std::vector<int> SMatrix::get_H_p_vectors(int p) {
    // TODO optimize < O(n) или сделать мапку с заранее сгенерированными векторами 
    std::vector<int> vecs;
    for (int i = 0; i < power_2_n(2*N); i++) {
        if (get_number_of_ones(i) == p) {
            vecs.push_back(i);
        }
    }

    return vecs;
}

int get_ones_positions(int n, int *pos1, int *pos2) {
    int nn = n;
    char was_one = 0;
    char was_two = 0;
    for (int i = 0; i < sizeof(int); i++) {
        int bit = nn & 1;
        if (bit) {
            if (!was_one) {
                *pos1 = i;
                was_one = 1;
            } else if (was_one && !was_two) {
                *pos2 = i;
                was_two = 1;
            } else {
                return -1;
            }
        }
        nn >>= 1;
    }

    return 0;
}

double SMatrix::get_interaction(int vec1, int vec2) {
    int x = vec1^vec2;
    int pos1, pos2;
    int res = get_ones_positions(x, &pos1, &pos2);
    if (res == -1) { // больше двух единиц
       return 0.0;

    } else if ((pos2 - pos1) > 2) { // индексы единиц слишком далеко друг от друга
        return 0.0;

    } else if ((pos2 - pos1) == 2) { 
        if (pos1 % 2 == 0) { // обмен фотонами между полостями
            return a;
        } else {
            return 0.0;
        }

    } else if ((pos2 - pos1) == 1) { 
        if ((pos2 / 2) == (pos1 / 2)) { // атом-фотон в одной полости
            return b;
        } else {
            return 0.0;
        }
    }

    return 0.0;
}


double **SMatrix::gen_H_p(int p, int sz) {
    double **H_p = new double*[sz];
    for (int i = 0; i < sz; i++) {
        H_p[i] = new double[sz];
    }

    // fill zeros
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            H_p[i][j] = 0.0;
        }
    }
    
    std::vector<int> H_p_vectors = get_H_p_vectors(p);

    for (int i = 0; i < sz; i++) {
        int vec1 = H_p_vectors[i];
        int n_of_fotons = get_number_of_ones_even(vec1);
        int n_of_atoms  = get_number_of_ones_odd(vec1);
        H_p[i][i] = n_of_atoms * w_a + n_of_fotons * w_c;

        for (int j = i+1; j < sz; j++) {
            int vec2 = H_p_vectors[j];
            H_p[i][j] = get_interaction(vec1, vec2);
        }
    }

    // отображаем относительно диагонали
    for (int i = 0; i < sz; i++) {
        for (int j = i + 1; j < sz; j++) {
            H_p[j][i] = H_p[i][j];
        }
    }

    return H_p;
}

double SMatrix::get_H_p_i_j(int p, int i, int j, std::map<int, double **> H_generated, std::map<int, int> H_sizes) {
    if (H_generated.find(p) == H_generated.end()) {
        H_generated[p] = gen_H_p(p, H_sizes[p]);
    }
    
    return H_generated[p][i][j];
}

void SMatrix::fill_RWA(std::vector<int> H_idxs, std::map<int, int> H_sizes) {
    barrier();

	std::map<int, double**> H_generated;

	int global_i;
    int global_j;

    int p;
    int H_i, H_j;

    for (int my_i = 0; my_i < myProcRows; my_i++) {
        global_i = myProcRowsOffset + my_i;

        for (int my_j = 0; my_j < myProcCols; my_j++) {
            global_j = myProcColsOffset + my_j;

            p = get_H_idx(global_i, global_j, H_idxs, H_sizes, &H_i, &H_j);
            if (p != -1) {
                // set my_H[my_i][my_i] = get_H_p_i_j(p, H_i, H_j);
				setij(my_i, my_j, get_H_p_i_j(p, H_i, H_j, H_generated, H_sizes));
            } else {
                // set my_H[my_i][my_i] = 0.0;
				setij(my_i, my_j, 0.0);
            }
        }
    }

	for (auto pair: H_generated) {
        double **to_delete = pair.second;
        for (int i = 0; i < H_sizes[pair.first]; i++) {
            delete[] to_delete[i];
        }

        delete[] to_delete;
    }


	barrier();
}