#include <iostream> 
#include <iomanip>
#include "smatrix.h"

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

/// SVD demo
void printSVD(SMatrix &x, bool isRoot) {
	int SSize;
	SMatrix *U, *VT;
	SReal *S = x.calculateSVD(&SSize, &U, &VT);

	if (isRoot) {
		cout << endl;
		cout << "Producing SVD:" << endl;
		cout << " =" << endl;
	}
	cout << *U;
	if (isRoot) {
		cout << " *" << endl;
		cout << "  diag([";
		for (int i = 0; i < SSize; i++) {
			if (i > 0) {
				cout << " ";
			}
			cout << setprecision(8) << S[i];
		}
		cout << "])" << endl;
		cout << " *" << endl;
	}
	cout << *VT;
	
	delete U, VT, SSize; 
	delete[] S;
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	bool isRoot = SMatrix::init(0, 0, 5);

	/*
#ifdef USE_COMPLEX
	{
		SMatrix x(10, 10);
		x.populate(&filler0);
		cout << x;
		printSVD(x, isRoot);
	}
#endif
	*/
	{
		SMatrix x(4, 5);
		x.setIdentity();
		x.populate(&fillerTest);
		cout << x;
		printSVD(x, isRoot);
	}
	if (isRoot) {
		cout << endl << endl << endl;
	}
	{
		SMatrix x(8, 8, 3, 1);
		x.dumpInfo();

		x = 0;
		SMatrix y = x;
		x.fill(12.1);
		
		cout << x;
		if (isRoot) {
			cout << endl;
		}

		x.setIdentity();
		cout << x;
		if (isRoot) {
			cout << endl;
		}
		
		x.populate(&filler);
		cout << x;
		if (isRoot) {
			cout << endl;
		}

		SType data[6] = {1, 2, 3, 4, 5, 6};
		y.barrier();
		y.set(data, 2, 3, 2, 3);
		y.set(data, 4, 2, 3, 2);
		y.barrier();
		cout << y;
		if (isRoot) {
			cout << endl;
		}
		{
			bool equals = (x == y);
			if (isRoot) {
				cout << (equals ? "Equal, BAD" : "Not equal, OK") << endl;
			}
		}
		x = y;
		cout << x;
		{
			bool equals = (x == y);
			if (isRoot) {
				cout << (equals ? "Equal, OK" : "Not equal, BAD") << endl;
			}
		}
		if (isRoot) {
			cout << endl;
		}

		SMatrix z(8, 8);
		z.populate(filler);
		cout << z;
		printSVD(z, isRoot);
	}
	if (isRoot) {
		cout << endl << endl << endl;
	}
	{
		SType data[] = {
			8.79,	9.93,	9.83,	5.45,	3.16,
			6.11,	6.91,	5.04,	-0.27,	7.98,
			-9.15,	-7.93,	4.86,	4.85,	3.01,
			9.57,	1.64,	8.83,	0.74,	5.80,
			-3.49,	4.02,	9.80,	10.00,	4.27,
			9.84,	0.15,	-8.99,	-6.02,	-5.31
			};
		SMatrix x(6, 5);
		x.set(data);
		cout << x;
		printSVD(x, isRoot);
	}

	SMatrix::exit();
	return 0;
}