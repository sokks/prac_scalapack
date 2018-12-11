#include "smatrix.h"

SMatrix::SMatrix() {
    
}

void trivial_evolution(int size, SMatrix& roStart, SMatrix& H, double dT, int n_steps) {
    SMatrix ro = roStart;
    SMatrix U = SMatrix::exp(ro * (complex_d(-1, 0) / complex_d(0,1) * complex_d(dT, 0)));

    for (int step = 0; step < n_steps; step++) {
        ro = SMatrix::conjugate(U) * ro * U;
        ro.print_diag();
    }
}