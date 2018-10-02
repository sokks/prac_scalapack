#pragma once
#include "../common.h"

complex_d filler(int i, int j) {
    return complex_d(i, j);
}

complex_d diag_filler(int i) {
    return complex_d(i, i+1.5);
}
