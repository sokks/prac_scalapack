#include "mwriter.h"
#include "filler.hpp"

complex_d test_filler(int i, int j)
{
    return complex_d(i, j);
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        std::cout << "usage: mwriter_test <file txt> <file diag txt> <file raw> <file diag raw>\n";
        return 0;
    }

    char *f_in = argv[1];
    char *f_in_diag = argv[2];
    char *f_out = argv[3];
    char *f_out_diag = argv[4];

    // convert full matrix txt to internal binary
    matrix_d m2;
    m2.read_from_txt(f_in);
    m2.write_to(f_out);
    m2 = matrix_d();
    m2.read_from(f_out);
    m2.print();

    // convert diagonal only txt to internal binary
    matrix_d m3;
    m3.read_diag_from_txt(f_in_diag);
    m3.write_diag_to(f_out_diag);
    m3 = matrix_d();
    m3.read_diag_from(f_out);
    m3.print_diag();

    // use custom filler for full matrix
    matrix_d m5(5);
    m5.fill(&filler);
    m5.print();

    // use custom diagonal filler
    matrix_d m6(5);
    m6.fill_diag(&diag_filler);
    m6.print();

    return 0;
}
