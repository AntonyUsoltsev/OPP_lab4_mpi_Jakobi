#include <iostream>
#include <vector>
#include <cmath>

constexpr double a = 1e5;
constexpr double eps = 1e-8;
constexpr double Dx = 2, Dy = 2, Dz = 2;
constexpr int N_x = 11, N_y = 11, N_z = 5;
constexpr double x_0 = -1, y_0 = -1, z_0 = -1;

double h_value(int D, int N) {
    return (double) D / (N - 1);
}

const double h_x = Dx / (N_x - 1);
const double h_y = Dy / (N_y - 1);
const double h_z = Dz / (N_z - 1);
//const double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + 2 / (h_z * h_z) + a);
const double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + a);
//const double coef = 1 / (2 / (h_x * h_x) + a);

class Grid_2D {
public:
    std::vector<double> grid;
    std::vector<double> next_phi;
    bool stop_flag = true;
    size_t grid_size = N_x * N_y;

    Grid_2D() {
        grid.resize(grid_size);
        next_phi.resize(grid_size);
    }

    void fill_grid() {
        for (int i = 0; i < grid_size; ++i) {
            grid[i] = 0;
        }
    }

    void count_next_phi(int i, int j) {

        double up_x_neughbor = (j * N_x + i + N_x) > (grid_size - 1) ? 1 : grid[j * N_x + i + N_x];
        double down_x_neughbor = (j * N_x + i - N_x) < 0 ? 1 : grid[j * N_x + i - N_x];
        double left_y_neighbor = ((j * N_x + i) % N_x == 0) ? 1 : grid[i - 1];
        double right_y_neighbor = ((j * N_x + i + 1) % N_x == 0) ? 1 : grid[i + 1];
        double coefc = coef;
        double v1 =
                (up_x_neughbor + down_x_neughbor) / (h_x * h_x) + (left_y_neighbor + right_y_neighbor) / (h_y * h_y);
        double phivar = phi(i, j, 0);
        double next_var = coefc * (v1 - 6 + a * phivar);
        next_phi[j * N_x + i] = next_var;

    }


    double phi(int i, int j, int k) {
        return (x_0 + i * h_x) * (x_0 + i * h_x) + (y_0 + j * h_y) * (y_0 + j * h_y);
    }

    void iteration() {
        for (int i = 0; i < N_x; ++i) {
            for (int j = 0; j < N_y; ++j) {
                count_next_phi(i, j);
            }
        }
        double maxdif = 0;
        for (int i = 0; i < grid_size; ++i) {
            double curdif = fabs(grid[i] - next_phi[i]);
            if (curdif > maxdif) {
                maxdif = curdif;
            }
            grid[i] = next_phi[i];
        }
        if (maxdif < eps) {
            stop_flag = false;
        }
    }

    void print_grid() {
        for (int i = 0; i < N_y; ++i) {
            for (int j = 0; j < N_x; ++j) {
                printf("%.2lf ", grid[i * N_x + j]);
            }
            std::cout << std::endl;
        }
    }

};


int main() {
    Grid_2D grid2D;
    grid2D.fill_grid();
    grid2D.print_grid();
    // std::cout << h_x << " " << N_x << " " << Dx << " " << coef<<" "<<(2 / (h_x * h_x) + a) ;
//    for (int i = 0; i < 10; i++) {
//        std::cout << "Iteration " << i+1 << '\n';
//        grid1D.iteration();
//        grid1D.print_grid();
//    }
    // std::cout << "\n\n";
    //  grid2D grid2D;
    //grid2D.fill_grid();
    int iter = 1;
    while (grid2D.stop_flag) {
        std::cout << "Iteration " << iter << '\n';
        grid2D.iteration();
        grid2D.print_grid();
        iter++;
    }
    std::cout << std::endl;
    for (int i = -10; i <= 10; i += 2) {
        for (int j = -10; j <= 10; j += 2) {
            double c = (double) j / 10;
            double b = (double) i / 10;
           // std::cout << c * c + b * b << " ";
            printf("%.2lf ",  c * c + b * b);
        }
        std::cout << std::endl;
    }


    return 0;
}
