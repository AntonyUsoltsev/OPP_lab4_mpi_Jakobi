#include <iostream>
#include <vector>
#include <cmath>

constexpr double a = 1e5;
constexpr double eps = 1e-8;
constexpr double Dx = 2, Dy = 2, Dz = 2;
constexpr int N_x = 5, N_y = 5, N_z = 5;
constexpr double x_0 = -1, y_0 = -1, z_0 = -1;

double h_value(int D, int N) {
    return (double) D / (N - 1);
}

const double h_x = Dx / (N_x - 1);
const double h_y = Dy / (N_y - 1);
const double h_z = Dz / (N_z - 1);
const double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + 2 / (h_z * h_z) + a);
//const double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + a);
//const double coef = 1 / (2 / (h_x * h_x) + a);

class Grid_3D {
public:
    std::vector<double> grid;
    std::vector<double> next_phi;
    bool stop_flag = true;
    size_t grid_size = N_x * N_y * N_z;

    Grid_3D() {
        grid.resize(grid_size);
        next_phi.resize(grid_size);
    }

    void fill_grid() {
        for (int i = 0; i < grid_size; ++i) {
            grid[i] = 0;
        }
    }

    void count_next_phi(int i, int j, int k) {
//        double up_y_neughbor = (cur_pos + N_x) > (grid_size - 1) ? 1 : grid[cur_pos + N_x];
//        double down_y_neughbor = (cur_pos - N_x) < 0 ? 1 : grid[cur_pos - N_x];
//
//        double left_x_neighbor = ((cur_pos) % N_x == 0) ? 1 : grid[i - 1];
//        double right_x_neighbor = ((cur_pos + 1) % N_x == 0) ? 1 : grid[i + 1];


        int cur_pos = k * N_x * N_y + j * N_x + i;

        double left_x_neighbor = (cur_pos % N_x == 0) ? 1 : grid[i - 1];
        double right_x_neighbor = ((cur_pos + 1) % N_x == 0) ? 1 : grid[i + 1];

        double backward_y_neighbor = (cur_pos - k * N_x * N_y + N_x) > (N_x * N_y - 1) ? 1 : grid[cur_pos + N_x];
        double forward_y_neighbor = (cur_pos - k * N_x * N_y - N_x) < 0 ? 1 : grid[cur_pos - N_x];

        double up_z_neighbor = (cur_pos + N_x * N_y) > (grid_size - 1) ? 1 : grid[cur_pos];;
        double down_z_neighbor = (cur_pos - N_x * N_y) < 0 ? 1 : grid[cur_pos];
//TODO:change N when parallel
        double coefc = coef;
        double v1 = (left_x_neighbor + right_x_neighbor) / (h_x * h_x) +
                    (backward_y_neighbor + forward_y_neighbor) / (h_y * h_y) +
                    (up_z_neighbor + down_z_neighbor) / (h_z * h_z);
        double phivar = phi(i, j, k);
        double next_var = coefc * (v1 - 6 + a * phivar);
        next_phi[cur_pos] = next_var;

    }

    double phi(int i, int j, int k) {
        return (x_0 + i * h_x) * (x_0 + i * h_x) + (y_0 + j * h_y) * (y_0 + j * h_y) +
               (z_0 + k * h_z) * (z_0 + k * h_z);
    }

    void iteration() {
        for (int i = 0; i < N_x; ++i) {
            for (int j = 0; j < N_y; ++j) {
                for (int k = 0; k < N_z; ++k) {
                    count_next_phi(i, j, k);
                }
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
        for (int k = 0; k < N_z; ++k) {
            std::cout << "Slice z =" << k<<std::endl;
            for (int j = 0; j < N_y; ++j) {
                for (int i = 0; i < N_x; ++i) {
                    printf("%.2lf ", grid[k * N_x * N_y + j * N_x + i]);
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

};


int main() {
    Grid_3D grid3D;
    grid3D.fill_grid();
    grid3D.print_grid();
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
    while (grid3D.stop_flag) {
        std::cout << "Iteration " << iter << '\n';
        grid3D.iteration();
        grid3D.print_grid();
        iter++;
    }
    std::cout << std::endl;

    for (int k = -10; k <= 10; k += 5) {
        std::cout << "Slice z =" << k << std::endl;
        for (int j = -10; j <= 10; j += 5) {
            for (int i = -10; i <= 10; i += 5) {
                double c = (double) j / 10;
                double b = (double) i / 10;
                double d = (double) k / 10;
                // std::cout << c * c + b * b << " ";
                printf("%.2lf ", c * c + b * b +d*d);
            }std::cout << std::endl;
        }
        std::cout << std::endl;
    }



    return 0;
}
