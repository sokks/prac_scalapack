#include "smatrix.h"

// usage: ./evol <size> <ro type> <ro file | номер базисного состояния> <dT> <H file> <n of steps>
// ro type: 1 - full matrix, 2 - diag, 3 - номер базисного состояния
int main(int argc, char **argv) {
    if (argc < 7) {
        std::cout << "usage: ./evol <size> <ro type> <ro file | номер базисного состояния> <dT> <H file> <n of steps>\n";
    }

    int    size    = std::atoi(argv[1]);
    char   ro_type = argv[2][0];
    char*  ro_file = argv[3];
    double dT      = std::stod(argv[4]);
    char*  H_file  = argv[5];
    int    n_steps = std::atoi(argv[6]);

    // todo
    SMatrix H;
    SMatrix Ro;

    if (ro_type == '1') {
        // todo 
        // SMatrix Ro(?, ro_file) -- распределенное чтение матрицы
    } else if (ro_type == '2') { 
        // todo 
        // SMatrix Ro(?, diag_file) -- распределенное чтение диагонали
    } else if (ro_type == '3') {
        // todo 
        // SMatrix Ro(2^size(?))
        // i = std::atoi(ro_file)
        // Ro.set(i,j, complex_d(1, 0))
    }

    trivial_evolution(size, Ro, H, dT, n_steps);

    return 0;
}