#include "../common.h"
#include <iostream>
#include <fstream>


class matrix_d {
    int N;
    double *data;

public:
    matrix_d(int N = 0);
    matrix_d(const matrix_d&);

    ~matrix_d();

    matrix_d& operator=(const matrix_d&);

    void fill(complex_d (*filler)(int, int));
    void fill_diag(complex_d (*diag_filler)(int));

    void read_from(const char *filename);
    void read_from_txt(const char *filename);
    void read_diag_from(const char *filename);
    void read_diag_from_txt(const char *filename);

    void write_to(const char *filename);
    void write_diag_to(const char *filename);

    void print();
    void print_diag();
};