#include <iostream>
#include <map>
#include <vector>
#include "smatrix.cpp"

int N = 0;

double w_c = 1.0;
double w_a = 10.0;
double a = 1.0;
double b = 100.0;

int fac(int n) {
    int res = 1;
    for (int i = 1; i <= n; i++) {
        res *= i;
    }

    return res;
}

int C_n_k(int n, int k) {
    return fac(n) / (fac(k) * fac(n-k));
}

int get_H_size(int H_idx, int N) {
    return C_n_k(2*N, H_idx);
}

void print_H(double **H, int sz) {
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            std::cout << H[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

        N     = std::atoi(argv[1]);
    int k     = std::atoi(argv[2]);
    int E_min = std::atoi(argv[3]);
    int E_max = std::atoi(argv[4]);

    double a   = std::atof(argv[5]);
    double b   = std::atof(argv[6]);
    double w_a = std::atof(argv[7]);
    double w_c = std::atof(argv[8]);

    std::vector<int> H_idxs;

    for (int s = 0; s <=k; s++) {
        for (int j = std::max(E_min - s, 0); j <= std::min(E_max - s, 2*N); j++) {
            H_idxs.push_back(j);
            std::cout << j << " ";
        }
    }

    std::cout << std::endl;

    std::map<int, int> H_sizes;
    for (int idx = 0; idx <= 2*N; idx++) {
        H_sizes[idx] = get_H_size(idx, N);
        std::cout << idx << " " << H_sizes[idx] << std::endl;
    }

    int H_size_full = 0;
    for (auto idx: H_idxs) {
        H_size_full += H_sizes[idx];
    }

    std::cout << H_size_full << std::endl;

    bool isRoot = SMatrix::init(0, 0, 2);
    int count, *buf;


    

    SMatrix x(H_size_full, H_size_full);
    x.N = N;
    x.w_a = w_a;
    x.w_c = w_c;
    x.a = a;
    x.b = b;

    x.fill(0);
    x.fill_RWA(H_idxs, H_sizes);

    // cout<<x;

    x.Print();

    // char *fileout = "data/gened_H.dat";
    // x.writef(fileout);

    SMatrix::exit();
    return 0;
}