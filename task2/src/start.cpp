#include <iostream>
#include <map>
#include <vector>

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
    
    // int up = 1;
    // int down = 1;

    // if (k == 0) {
    //     return 1;
    // }

    // for (int i = n-k+1; i <= n; i++) {
    //     up *= i;
    // }

    // for (int i = 1; i <= n-k; i++) {
    //     down *= i;
    // }

    // std::cout << up << " " << down << std::endl;
    // return up / down;
}

int get_H_size(int H_idx, int N) {
    return C_n_k(2*N, H_idx);
}

int get_H_idx(int global_i, int global_j, std::vector<int> H_idxs, std::map<int, int> H_sizes, int *H_i, int *H_j) {
    int idx_i = -1;
    int idx_j = -1;
    
    int prev_sum_i = 0;
    for (int i = 0; i < H_idxs.size(); i++) {
        int sz = H_sizes[H_idxs[i]];
        if (global_i < prev_sum_i + sz) {
            idx_i = H_idxs[i];
            break;
        }
        prev_sum_i += sz;
    }

    int prev_sum_j = 0;
    for (int j = 0; j < H_idxs.size(); j++) {
        int sz = H_sizes[H_idxs[j]];
        if (global_j < prev_sum_j + sz) {
            idx_j = H_idxs[j];
            break;
        }
        prev_sum_j += sz;
    }

    if ( (idx_i != -1) && (idx_j != -1) && (idx_i == idx_j) ) {
        *H_i = global_i - prev_sum_i;
        *H_j = global_j - prev_sum_j;
        return idx_i;
    } 

    return -1;
}



int get_number_of_ones(int n) {
    unsigned int count = 0;
    int n_start = n;
    for (; n; n >>= 1) {
        count += n & 1;
    }

    return count;
}

int get_number_of_ones_odd(int n) {
    unsigned int count = 0;
    for (n = n >> 1; n; n >>= 2) {
        count += n & 1;
    }
    return count;
}

// количество единиц на четных местах (если считать с конца с 0) 
// н-р: [00...00101 -> 2] [00...00010 -> 0]
int get_number_of_ones_even(int n) {
    unsigned int count = 0;
    for (; n; n >>= 2) {
        count += n & 1;
    }
    return count;
}

int power_2_n(int n) {
    int res = 1;
    for (int i = 0; i < n; i++) {
        res <<= 1;
    }
    return res;
}

std::vector<int> get_H_p_vectors(int p) {
    // TODO optimize < O(n)
    std::vector<int> vecs;
    for (int i = 0; i < power_2_n(2*N); i++) {
        if (get_number_of_ones(i) == p) {
            vecs.push_back(i);
        }
    }

    return vecs;
}

int get_ones_positions(int n, int *pos1, int *pos2) {
    int nn = n;
    char was_one = 0;
    char was_two = 0;
    for (int i = 0; i < sizeof(int); i++) {
        int bit = nn & 1;
        if (bit) {
            if (!was_one) {
                *pos1 = i;
                was_one = 1;
            } else if (was_one && !was_two) {
                *pos2 = i;
                was_two = 1;
            } else {
                return -1;
            }
        }
        nn >>= 1;
    }

    return 0;
}

double get_interaction(int vec1, int vec2) {
    int x = vec1^vec2;
    int pos1, pos2;
    int res = get_ones_positions(x, &pos1, &pos2);
    if (res == -1) { // больше двух единиц
       return 0.0;

    } else if ((pos2 - pos1) > 2) { // индексы единиц слишком далеко друг от друга
        return 0.0;

    } else if ((pos2 - pos1) == 2) {
        if (pos1 % 2 == 0) {
            return a;
        } else {
            return 0.0;
        }

    } else if ((pos2 - pos1) == 1) { 
        if ((pos2 / 2) == (pos1 / 2)) {
            return b;
        } else {
            return 0.0;
        }
    }

    return 0.0;
}

double **gen_H_p(int p, int sz) {
    double **H_p = new double*[sz];
    for (int i = 0; i < sz; i++) {
        H_p[i] = new double[sz];
    }

    // fill zeros
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            H_p[i][j] = 0.0;
        }
    }
    
    std::vector<int> H_p_vectors = get_H_p_vectors(p);

    for (int i = 0; i < sz; i++) {
        int vec1 = H_p_vectors[i];
        int n_of_fotons = get_number_of_ones_even(vec1);
        int n_of_atoms  = get_number_of_ones_odd(vec1);
        H_p[i][i] = n_of_atoms * w_a + n_of_fotons * w_c;

        for (int j = i+1; j < sz; j++) {
            int vec2 = H_p_vectors[j];
            H_p[i][j] = get_interaction(vec1, vec2);
        }
    }

    // отображаем относительно диагонали
    for (int i = 0; i < sz; i++) {
        for (int j = i + 1; j < sz; j++) {
            H_p[j][i] = H_p[i][j];
        }
    }

    return H_p;
}

double get_H_p_i_j(int p, int i, int j, std::map<int, double **> H_generated, std::map<int, int> H_sizes) {
    if (H_generated.find(p) == H_generated.end()) {
        H_generated[p] = gen_H_p(p, H_sizes[p]);
    }
    
    return H_generated[p][i][j];
}

void print_H(double **H, int sz) {
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            std::cout << H[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
        N     = 2;
    int E_min = 0;
    int E_max = 4;
    int k     = 1;

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

    int proc_n_i = 0;
    int proc_n_j = 0;

    int proc_sz_i = 10;
    int proc_sz_j = 10;

    int proc_start_i = proc_n_i * proc_sz_i;
    int proc_start_j = proc_n_j * proc_sz_j;

    // int my_i = 2;
    // int my_j = 3;

    int global_i;
    int global_j;

    int p;
    int H_i, H_j;

    for (int my_i = 0; my_i < proc_sz_i; my_i++) {
        global_i = proc_start_i + my_i;

        for (int my_j = 0; my_j < proc_sz_j; my_j++) {
            global_j = proc_start_j + my_j;

            p = get_H_idx(global_i, global_j, H_idxs, H_sizes, &H_i, &H_j);
            if (p != -1) {
                // set my_H[my_i][my_i] = get_H_p_i_j(p, H_i, H_j);
            } else {
                // set my_H[my_i][my_i] = 0.0;
            }
        }
    }

    std::map<int, double**> H_generated;

    for (auto idx: H_idxs) {
        if (H_generated.find(idx) == H_generated.end()) {
            std::cout << "generating H with idx " << idx << std::endl;
            H_generated[idx] = gen_H_p(idx, H_sizes[idx]);
        }

        std::cout << "H_idx: " << idx << std::endl;
        print_H(H_generated[idx], H_sizes[idx]);
    }


    for (auto pair: H_generated) {
        double **to_delete = pair.second;
        for (int i = 0; i < H_sizes[pair.first]; i++) {
            delete[] to_delete[i];
        }

        delete[] to_delete;
    }

    return 0;
}