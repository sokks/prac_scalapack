#include "mwriter.h"
#include "filler.hpp"

complex_d test_filler(int i, int j)
{
    return complex_d(i, j);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "usage: mwriter_test <size> <filename>\n";
        return 0;
    }

    int size = std::atoi(argv[1]);
    char *filename = argv[2];

    matrix_d m(size, size);
    m.fill(&filler);
    m.write_to(filename);

    matrix_d m_check;
    m_check.read_from(filename);
    m_check.print();
    return 0;
}
