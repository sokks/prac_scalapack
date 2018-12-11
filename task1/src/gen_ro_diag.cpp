#include "mwriter/mwriter.h"
#include "common.h"

#include <ctime>

int seed = 1;
complex_d *diag = nullptr;

void gen_diag(int N) {
    diag = new complex_d[N];
    std::srand(seed);
    double sum = 0;
    for (int i = 0; i < N; i++) {
        double lambda = double(std::rand());
        diag[i] = complex_d( lambda * lambda, 0);
        sum += diag[i].real();
    }
    for (int i = 0; i < N; i++) {
        diag[i] /= sum;
    }
}

complex_d diag_filler(int i) {
    return diag[i];
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "usage: gen_ro_diag <filename> <N> [<seed>]\n";
        return 0;
    }

    char *filename = argv[1];
    int N = std::atoi(argv[2]);
    if (argc > 3) {
        seed = std::atoi(argv[3]);
    } else {
        seed = std::time(0);
    }

    gen_diag(N);
    matrix_d m_diag(N);
    m_diag.fill_diag(&diag_filler);
    m_diag.print_diag();
    m_diag.write_diag_to(filename);

    return 0;
}
