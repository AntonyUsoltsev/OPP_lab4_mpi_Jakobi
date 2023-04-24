#include <iostream>
#include <vector>
#include <cmath>

constexpr double a = 1e5;
constexpr double eps = 1e-8;
constexpr double Dx = 2, Dy = 2, Dz = 2;
constexpr int N_x = 11, N_y = 11, N_z = 11;
constexpr double x_0 = -1, y_0 = -1, z_0 = -1;

constexpr double h_x = Dx / (N_x - 1);
constexpr double h_y = Dy / (N_y - 1);
constexpr double h_z = Dz / (N_z - 1);
constexpr double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + 2 / (h_z * h_z) + a);


class Grid_3D {
public:
    std::vector<double> grid;
    std::vector<double> next_phi;
    bool stop_flag = true;

    int x_size = N_x + 2, y_size = N_y, z_size = N_z;

    size_t grid_size = x_size * y_size * z_size;


    double phi(int i, int j, int k) {
        return (x_0 + i * h_x) * (x_0 + i * h_x) + (y_0 + j * h_y) * (y_0 + j * h_y) +
               (z_0 + k * h_z) * (z_0 + k * h_z);
    }

    int grid_ind_from_coords(int x, int y, int z) {
        return y_size * z_size * x + y + y_size * z;
    }

    Grid_3D(int x_size, int y_size, int z_size)  {
        this->x_size = x_size;
        this->y_size = y_size;
        this->z_size = z_size;
        grid.resize(grid_size);
        next_phi.resize(grid_size);
    }

    void fill_grid() {
        std::fill(grid.begin(), grid.end(), 0);
        std::fill(next_phi.begin(), next_phi.end(), 0);
        for (int i = 0; i < z_size; i++) {
            for (int j = 0; j < y_size; j++) {
                grid[grid_ind_from_coords(0, j, i)] = phi(-1, j, i);
                grid[grid_ind_from_coords(x_size - 1, j, i)] = phi(x_size - 2, j, i);
            }
        }
    }



    void count_next_phi(int i, int j, int k) {

        int cur_pos = grid_ind_from_coords(i, j, k);

        double left_x_neighbor = grid[grid_ind_from_coords(i - 1, j, k)];
        double right_x_neighbor = grid[grid_ind_from_coords(i + 1, j, k)];

        double backw_y_neighbor = (j >= (y_size - 1)) ? phi(i - 1, j + 1, k) : grid[grid_ind_from_coords(i, j + 1, k)];
        double forw_y_neighbor = (j <= 0) ? phi(i - 1, -1, k) : grid[grid_ind_from_coords(i, j - 1, k)];

        double up_z_neighbor = (k >= (z_size - 1)) ? phi(i - 1, j, k + 1) : grid[grid_ind_from_coords(i, j, k + 1)];
        double down_z_neighbor = (k <= 0) ? phi(i - 1, j, -1) : grid[grid_ind_from_coords(i, j, k - 1)];

////TODO:change N when parallel
        double v1 = (left_x_neighbor + right_x_neighbor) / (h_x * h_x) +
                    (backw_y_neighbor + forw_y_neighbor) / (h_y * h_y) +
                    (up_z_neighbor + down_z_neighbor) / (h_z * h_z);
        double phivar = phi(i - 1, j, k);
        double next_var = coef * (v1 - 6 + a * phivar);
        next_phi[cur_pos] = next_var;

    }


    void iteration() {
        for (int i = 1; i < x_size - 1; ++i) {
            for (int k = 0; k < z_size; ++k) {
                for (int j = 0; j < y_size; ++j) {
                    count_next_phi(i, j, k);
                }
            }
        }
        double maxdif = 0;
        for (int i = y_size * z_size; i < grid_size - y_size * z_size; ++i) {
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
        for (int k = 0; k < x_size; ++k) {
            std::cout << "Slice x = " << k << std::endl;
            for (int i = 0; i < z_size; i++) {
                for (int j = 0; j < y_size; j++) {
                    printf("%.2lf ", grid[grid_ind_from_coords(k, i, j)]);
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    void print_grid(FILE *file) {
        for (int k = 1; k < x_size - 1; ++k) {
            fprintf(file, "Slice x = %.2lf\n", (x_0 + (k - 1) * h_x));
            for (int i = 0; i < z_size; i++) {
                for (int j = 0; j < y_size; j++) {
                    fprintf(file, "%.2lf ", grid[grid_ind_from_coords(k, i, j)]);
                }
                fputc('\n', file);
            }
            fputc('\n', file);
        }
    }

};


int main() {
    Grid_3D grid3D(7,N_y,N_z);
   // grid3D.x_size =7;
    grid3D.fill_grid();
    int iter = 1;
    while (grid3D.stop_flag) {
        std::cout << "Iteration " << iter << '\n';
        grid3D.iteration();
        iter++;
    }
    FILE *file = fopen("./../result.txt", "w");
    grid3D.print_grid(file);
    fclose(file);


    std::cout << std::endl;

//    for (int k = -10; k <= 10; k += 5) {
//        std::cout << "Slice x =" << k << std::endl;
//        for (int j = -10; j <= 10; j += 5) {
//            for (int i = -10; i <= 10; i += 5) {
//                double c = (double) j / 10;
//                double b = (double) i / 10;
//                double d = (double) k / 10;
//                // std::cout << c * c + b * b << " ";
//                printf("%.2lf ", c * c + b * b + d * d);
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }


    return 0;
}

//Slice x = 0
//4.25 3.50 3.25 3.50 4.25
//3.50 2.75 2.50 2.75 3.50
//3.25 2.50 2.25 2.50 3.25
//3.50 2.75 2.50 2.75 3.50
//4.25 3.50 3.25 3.50 4.25