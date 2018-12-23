#include <iostream>
#include <map>
#include <vector>

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

double get_H_p_i_j(int p, int i, int j) {
    // TODO
    
    return 0.0;
}

int main() {
    int N     = 4;
    int E_min = 3;
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


    return 0;
}