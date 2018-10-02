#include "mwriter.h"


matrix_d::matrix_d(int N, int M)
{
    this->N = N;
    this->M = M;
    this->data = nullptr;

    if (N && M)
    {
        this->data = new double[N * M * 2];
        for (int i = 0; i < N * M * 2; i++)
        {
            this->data[i] = 0.0;
        }
    }
}

matrix_d::matrix_d(matrix_d &m)
{
    this->N = m.N;
    this->M = m.M;

    if (this->data != nullptr)
    {
        delete[] this->data;
        this->data = nullptr;
    }

    if (m.data != nullptr)
    {
        this->data = new double[this->N * this->M * 2];
        for (int i = 0; i < this->N * this->M * 2; i++)
        {
            this->data[i] = m.data[i];
        }
    }
}

matrix_d::~matrix_d()
{
    if (data != nullptr)
    {
        free(data);
    }
}

void matrix_d::fill(complex_d (*filler)(int, int))
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            complex_d val = filler(i, j);
            this->data[(i * M + j) * 2] = val.real();
            this->data[(i * M + j) * 2 + 1] = val.imag();
        }
    }
}

void matrix_d::write_to(const char *filename)
{
    std::ofstream out(filename, std::ios::binary);
    out.write((char *)(&(this->N)), sizeof(int));
    out.write((char *)(&(this->M)), sizeof(int));
    out.write((char *)(this->data), this->N * this->M * 2 * sizeof(double));
    out.close();
}

void matrix_d::read_from(const char *filename)
{
    std::ifstream in(filename, std::ios::binary);
    in.read((char *)(&(this->N)), sizeof(int));
    in.read((char *)(&(this->M)), sizeof(int));
    if (this->data != nullptr)
    {
        delete[] this->data;
    }
    this->data = new double[this->N * this->M];
    in.read((char *)(this->data), this->N * this->M * 2 * sizeof(double));
    in.close();
}

void matrix_d::print()
{
    std::cout.setf(std::ios::fixed);
    std::cout.precision(4);
    std::cout << this->N << " " << this->M << std::endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            printf("( %3.3f, %3.3f) ", this->data[(i * M + j) * 2], this->data[(i * M + j) * 2 + 1]);
        }
        std::cout << std::endl;
    }
}
