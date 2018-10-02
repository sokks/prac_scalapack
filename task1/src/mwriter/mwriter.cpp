#include "mwriter.h"

matrix_d::matrix_d(int N)
{
    this->N = N;
    this->data = nullptr;

    if (N)
    {
        this->data = new double[N * N * 2];
        for (int i = 0; i < N * N * 2; i++)
        {
            this->data[i] = 0.0;
        }
    }
}

matrix_d::matrix_d(const matrix_d &m)
{
    this->N = m.N;
    this->data = nullptr;

    if (m.data != nullptr)
    {
        this->data = new double[this->N * this->N * 2];
        for (int i = 0; i < this->N * this->N * 2; i++)
        {
            this->data[i] = m.data[i];
        }
    }
}

matrix_d &matrix_d::operator=(const matrix_d &m)
{
    this->N = m.N;

    if (this->data != nullptr)
    {
        delete[] this->data;
        this->data = nullptr;
    }

    if (m.data != nullptr)
    {
        this->data = new double[this->N * this->N * 2];
        for (int i = 0; i < this->N * this->N * 2; i++)
        {
            this->data[i] = m.data[i];
        }
    }

    return *this;
}

matrix_d::~matrix_d()
{
    if (data != nullptr)
    {
        delete[] data;
    }
}

void matrix_d::fill(complex_d (*filler)(int, int))
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            complex_d val = filler(i, j);
            this->data[(i * N + j) * 2] = val.real();
            this->data[(i * N + j) * 2 + 1] = val.imag();
        }
    }
}

void matrix_d::fill_diag(complex_d (*diag_filler)(int))
{
    for (int i = 0; i < N; i++)
    {
        complex_d val = diag_filler(i);
        this->data[(i * N + i) * 2] = val.real();
        this->data[(i * N + i) * 2 + 1] = val.imag();
    }
}

void matrix_d::write_to(const char *filename)
{
    std::ofstream out(filename, std::ios::binary);
    out.write((char *)(&(this->N)), sizeof(int));
    out.write((char *)(this->data), this->N * this->N * 2 * sizeof(double));
    out.close();
}

void matrix_d::write_diag_to(const char *filename)
{
    double* diag = new double[this->N * 2];
    for (int i = 0; i < this->N; i++) {
        diag[i*2] = this->data[(i*N+i)*2];
        diag[i*2+1] = this->data[(i*N+i)*2+1];
    }

    std::ofstream out(filename, std::ios::binary);
    out.write((char *)(&(this->N)), sizeof(int));
    out.write((char *)diag, this->N * 2 * sizeof(double));
    out.close();

    delete[] diag;
}

void matrix_d::read_from(const char *filename)
{
    std::ifstream in(filename, std::ios::binary);
    in.read((char *)(&(this->N)), sizeof(int));
    
    if (this->data != nullptr)
    {
        delete[] this->data;
    }
    this->data = new double[this->N * this->N * 2];

    in.read((char *)(this->data), this->N * this->N * 2 * sizeof(double));
    in.close();
}

void matrix_d::read_from_txt(const char *filename)
{
    std::ifstream in(filename);
    in >> (this->N);
    
    if (this->data != nullptr)
    {
        delete[] this->data;
    }
    this->data = new double[this->N * this->N * 2];

    for (int i = 0; i < this->N * this->N * 2; i++) {
        in >> this->data[i];
    }

    in.close();
}

void matrix_d::read_diag_from(const char *filename)
{
    std::ifstream in(filename, std::ios::binary);
    in.read((char *)(&(this->N)), sizeof(int));
    double* diag = new double[this->N * 2];
    in.read((char *)diag, this->N * this->N * 2 * sizeof(double));
    in.close();

    if (this->data != nullptr)
    {
        delete[] this->data;
    }
    this->data = new double[this->N * this->N * 2];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                this->data[(i*N+i)*2] = diag[i*2];
                this->data[(i*N+i)*2+1] = diag[i*2+1];
            } else {
                this->data[(i*N+i)*2] = 0.0;
                this->data[(i*N+i)*2+1] = 0.0;
            }
        }
    }
}

void matrix_d::read_diag_from_txt(const char *filename)
{
    std::ifstream in(filename);
    in >> this->N;

    double* diag = new double[this->N * 2];
    for (int i = 0; i < this->N*2; i++) {
        in >> diag[i];
    }

    in.close();

    if (this->data != nullptr)
    {
        delete[] this->data;
    }
    this->data = new double[this->N * this->N * 2];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                this->data[(i*N+i)*2] = diag[i*2];
                this->data[(i*N+i)*2+1] = diag[i*2+1];
            } else {
                this->data[(i*N+i)*2] = 0.0;
                this->data[(i*N+i)*2+1] = 0.0;
            }
        }
    }
}

void matrix_d::print()
{
    std::cout.setf(std::ios::fixed);
    std::cout << this->N << std::endl;
    for (int i = 0; i < this->N; i++)
    {
        for (int j = 0; j < this->N; j++)
        {
            printf("%3.4f %3.4f     ", this->data[(i * this->N + j) * 2], this->data[(i * this->N + j) * 2 + 1]);
        }
        std::cout << std::endl;
    }
}

void matrix_d::print_diag()
{
    std::cout.setf(std::ios::fixed);
    std::cout << this->N << std::endl;
    for (int i = 0; i < this->N; i++)
    {
        printf("%3.4f %3.4f     ", this->data[(i * this->N + i) * 2], this->data[(i * this->N + i) * 2 + 1]);
    }
    std:: cout << std::endl;
}
