#include "mwriter/mwriter.h"
#include "common.h"


int seed = 100;

int N = 0;
complex_d *full = nullptr;

void normalize() {
    complex_d full_sum(0,0);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            full_sum += full[i*N+j] * full[j*N+i];
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            full_sum += full[i*N+j] * full[j*N+i];
        }
    }

}
void gen_full(int N) {
    full = new complex_d[N*N];
    std::srand(seed);
    double sum = 0;
    for (int i = 0; i < N; i++) {
        full[i*N+i] = complex_d( double(std::rand()) / std::rand(), 0);
        sum += full[i*N+i].real();
        for (int j = i+1; j < N; j++) {
            full[i*N+j] = complex_d( double(std::rand()) / std::rand(), double(std::rand()) / std::rand());
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            full[i*N+j] = complex_d(full[j*N+i].real(), -full[j*N+i].imag());
        }
    }
}

complex_d full_filler(int i, int j) {
    return full[i*N+j];
    // if ( (abs(i-j) == 1) && (i>0) && (j>0) ) {
    //     return 1;
    // } 
    // if ((i == 2) && (j == 3)) {
    //     return 1;
    // }
    // if ((i == 3) && (j == 2)) {
    //     return 1;
    // }
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "usage: gen_H <filename> <N> [<seed>]\n";
        return 0;
    }

    char *filename = argv[1];
    N = std::atoi(argv[2]);

    // if (argc > 3) {
        seed = std::atoi(argv[3]);
    // } else {
        // seed = std::time(0);
    // }

    gen_full(N);
    matrix_d m_full(N);
    m_full.fill(&full_filler);
    m_full.print();
    m_full.write_to(filename);

    return 0;
}
