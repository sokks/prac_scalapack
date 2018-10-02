#include "../common.h"
#include <iostream>
#include <fstream>


class matrix_d {
    int N, M;
    double *data;

public:
    matrix_d(int N = 0, int M = 0);
    matrix_d(matrix_d&);

    ~matrix_d();

    matrix_d& operator=(matrix_d&);

    void fill(complex_d (*filler)(int, int));
    void write_to(const char *filename);
    void read_from(const char *filename);
    void print();
};