#include <iostream>
#include <vector>

constexpr int a = 10e5;
constexpr double eps = 10e-8;
constexpr int Dx = 2, Dy = 2, Dz = 2;
constexpr int Nx = 5, Ny = 5, Nz = 5;

double h_value(int D, int N) {
    return (double) D / (N - 1);
}

const double h_x = h_value(Dx, Nx);
const double h_y = h_value(Dy, Ny);
const double h_z = h_value(Dz, Nz);
const double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + 2 / (h_z * h_z) + a);

class Grid_3D {
public:
    std::vector<double> cur_phi;
    std::vector<double> next_phi;

    Grid_3D() {
        cur_phi.resize(Nx * Ny);
        next_phi.resize(Nx * Ny);
    }

    void fill_matrix() {
        for (int i = 0; i < Nx * Ny; ++i) {
            cur_phi[i] = 0;
        }
    }

    void count_next_phi(int i, int j) {
        double var1 = (j * Nx + i + Nx) > (Nx * Ny - 1) ? 1 : cur_phi[j * Nx + i + Nx];
        double var2 = (j * Nx + i - Nx) < 0 ? 1 : cur_phi[j * Nx + i - Nx];


    }


    int phi(int i, int j, int k) {
        return i * i + j * j + k * k;
    }

};


int main() {
    std::cout << a << std::endl;
    return 0;
}
