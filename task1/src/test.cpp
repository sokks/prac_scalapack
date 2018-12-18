#include <iostream>
#include <iomanip>
#include "matrix/smatrix.h"
#include <complex>

using namespace std;

SType filler(int i, int j) {
	i++;
	j++;
	return (i == j) ? i : (0.1*j+0.001*i);
}
SType fillerTest(int i, int j) {
	if (i == 0 && j == 0) {
		return 1;
	} else if (i == 0 && j == 4) {
		return 2;
	} else if (i == 1 && j == 2) {
		return 3;
	} else if (i == 3 && j == 1) {
		return 4;
	}
	return 0;
}
#ifdef USE_COMPLEX
SType filler0(int i, int j) {
	i++;
	j++;
	return SType(10 * i + j, 10 * j + i);
}
#endif

void printEIGENVAL(SMatrix &x, bool isRoot) {
	int SSize;
	SMatrix *Z;
	SReal *W = x.calculateEIGENVAL(&SSize, &Z);

	if (isRoot) {
		cout << endl;
		cout << "Producing Eigenvector:" << endl;
		cout << " =" << endl;
	}
	cout << *Z;
	if (isRoot) {
		cout << " *" << endl;
		cout << "  diag([";
		for (int i = 0; i < SSize; i++) {
			if (i > 0) {
				cout << " ";
			}
			cout << setprecision(8) << W[i];
		}
		cout << "])" << endl;
		cout << " *" << endl;
	}

	delete Z, SSize;
	delete[] W;
}

SType * diagonal_conversion(SReal *W, int SSize, double deltat, bool isRoot) {
	SType *W1 = new SType[SSize];

	for (int j = 0; j < SSize; j++) {
			W1[j] = exp(-W[j] * complex<SReal>(0,1) * deltat);
			// W1[j] = W[j];
		}
	return W1;
}


int main(int argc, char **argv) {
	
	MPI_Init(&argc, &argv);

	int SSize;

	/*
	 * Чтение гамильтониана и вычисление UdT
	 */

	char *file_H = argv[4];
	double deltat = atof(argv[3]);

	bool isRoot = SMatrix::init(0, 0, 2);
	int count, *buf;
	

	MPI_Status status;
	MPI_Offset filesize;
	MPI_File thefile, wrfile;
    
	MPI_File_open(MPI_COMM_WORLD, file_H, MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);
	MPI_File_set_view(thefile, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
	
	// read matrix sizes
	buf = (int *) malloc (2 * sizeof(int));
	SSize = buf[0]; // считаем матрицу квадратной
	MPI_File_read(thefile, buf, 2, MPI_INT, &status);
	// create matrix and read separately
	SMatrix x(buf[0], buf[1]);
	x.readf(thefile);

	/*SType data[] = {
		0, 2, -8,
		2, 0, 10,
		-8, 10, 5
		};*/
	//x.set(data);
	

	SMatrix *A,                                 // здесь будут собственные векторы по столбцам
			*D   = new SMatrix(SSize, SSize),   // здесь будет exp(-i * dT * eigen(H))
			*U   = new SMatrix(SSize, SSize),   // здесб будет UdT
			*tmp = new SMatrix(SSize, SSize);   // промежуточный результат

	SReal *W = x.calculateEIGENVAL(&SSize, &A);
	SType *W1 = diagonal_conversion(W, SSize, deltat, isRoot);
	*D = *A;
	(*D).DiagToMat(&SSize, W1);
	cout << *D << endl << endl; 
	cout << *A << endl << endl;

	(*A).mul((char *) "N", (char *) "N", D, &tmp);
	cout<< *tmp << endl << endl;

	SMatrix *UdT = new SMatrix(SSize, SSize);
	(*tmp).mul((char *) "N", (char *) "C", A, &UdT);
	cout<< *UdT << endl << endl;
	
	MPI_File_close(&thefile);
	
	
	/*
	 * Чтение матрицы плотности и вычисление
	 */
	SType wrbuf[SSize], wrbuf0[SSize];
	for (int i = 0; i < SSize; i++) {
		wrbuf[i] = complex<SReal>(0,0);
	}
	if (atoi(argv[2]) == 0)
	{
		MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile);
		MPI_File_set_view(thefile, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
		if (buf[1] != 0){
			x.readf(thefile);
			*A = x;
		}
		else {
			SMatrix v(buf[0], 1);
			v.readf(thefile);
			v.mul((char *) "C", (char *) "N", &v, &A);
		}
		
		
		MPI_File_read(thefile, buf, 2, MPI_INT, &status);
		MPI_File_close(&thefile);
	}
	else {
		x.fill(0);
		x.setij(atoi(argv[2])-1, atoi(argv[2])-1);
		*A = x;
	}
	
	// now A contains Ro
	cout << *A << endl << endl;
	int n = atoi(argv[5]);
	for (int j = 0; j < n; j++){
		(*UdT).mul((char *) "C", (char *) "N", A, &tmp);
		(*tmp).mul((char *) "N", (char *) "N", UdT, &A);
		(*A).wrbuf(wrbuf);
		MPI_Reduce(wrbuf, wrbuf0, SSize, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (isRoot) {
			for (int i = 0; i < SSize; i++) cout << wrbuf0[i];
			cout<<'\n';
		}
	}
	delete A, D, UdT, tmp;

	SMatrix::exit();
	return 0;
}
