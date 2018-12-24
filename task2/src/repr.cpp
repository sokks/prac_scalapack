#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
    char *filename = argv[2];
    int sz;

    std::ifstream in(filename, std::ios::binary);
    in.read((char *)(&sz), sizeof(int));
    std::cout << sz << std::endl;
    double *data = new double[sz];
    in.read((char *)(data), sz * sz * sizeof(double));
    in.close();

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            std::cout << data[i*sz+j] << " ";
        }
        std::cout << std::endl;
    }
}