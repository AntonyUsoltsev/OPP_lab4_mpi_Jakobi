#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

constexpr int RANK_ROOT = 0;

constexpr double a = 1e5;
constexpr double eps = 1e-8;
constexpr double Dx = 2.0f, Dy = 2.0f, Dz = 2.0f;
constexpr int N_x = 6, N_y = 6, N_z = 6;
constexpr double x_0 = -1.0f, y_0 = -1.0f, z_0 = -1.0f;

constexpr double h_x = Dx / (N_x - 1);
constexpr double h_y = Dy / (N_y - 1);
constexpr double h_z = Dz / (N_z - 1);
constexpr double coef = 1.0f / (2.0f / (h_x * h_x) + 2.0f / (h_y * h_y) + 2.0f / (h_z * h_z) + a);


class Grid_3D {
public:
    std::vector<double> grid;
    std::vector<double> next_phi;

    double maxdif, cur_maxdif;

    int grid_size, x_size, y_size, z_size;

    int comm_rank;

    double phi(int i, int j, int k) {
        int x_displ = i + comm_rank * (x_size - 2);
        return (x_0 + x_displ * h_x) * (x_0 + x_displ * h_x) + (y_0 + j * h_y) * (y_0 + j * h_y) +
               (z_0 + k * h_z) * (z_0 + k * h_z);
    }

    int grid_ind_from_coords(int x, int y, int z) const {
        return y_size * z_size * x + y + y_size * z;
    }

    Grid_3D(int x_size, int y_size, int z_size, int rank) {
        if (x_size > N_x or y_size > N_y or z_size > N_z) {
            throw std::runtime_error("Bad size");
        }
        this->x_size = x_size + 2;
        this->y_size = y_size;
        this->z_size = z_size;
        this->comm_rank = rank;
        maxdif = 100;
        printf("X size: %d, y size: %d, z size: %d, rank: %d\n", x_size, y_size, z_size, rank);
        std::cout.flush();
        grid_size = this->x_size * y_size * z_size;
        grid.resize(grid_size);
        next_phi.resize(grid_size);
    }

    void fill_grid(int size) {
        std::fill(grid.begin(), grid.end(), 0);
        std::fill(next_phi.begin(), next_phi.end(), 0);

        if (comm_rank == RANK_ROOT) {
            for (int i = 0; i < z_size; i++) {
                for (int j = 0; j < y_size; j++) {
                    grid[grid_ind_from_coords(0, j, i)] = phi(-1, j, i);
                }
            }
        }
        if (comm_rank == (size - 1)) {
            for (int i = 0; i < z_size; i++) {
                for (int j = 0; j < y_size; j++) {
                    grid[grid_ind_from_coords(x_size - 1, j, i)] = phi(x_size - 2, j, i);
                }
            }
        }
    }

    void iteration(int comm_size) {
        int layer_size = y_size * z_size;

        MPI_Request request_send_left, request_send_right, request_recv_left, request_recv_right;
        double *right_send_border = grid.data() + z_size * y_size * (x_size - 2);
        double *right_recv_border = grid.data() + z_size * y_size * (x_size - 1);
        if (comm_size > 1) {
            if (comm_rank == RANK_ROOT) {
                MPI_Isend(right_send_border, layer_size, MPI_DOUBLE, comm_rank + 1, 1, MPI_COMM_WORLD,&request_send_right);
            } else if (comm_rank == comm_size - 1) {
                MPI_Isend(grid.data() + z_size * y_size, layer_size, MPI_DOUBLE, comm_rank - 1, 1, MPI_COMM_WORLD,&request_send_left);
            } else {
                MPI_Isend(right_send_border, layer_size, MPI_DOUBLE, comm_rank + 1, 1, MPI_COMM_WORLD,&request_send_right);
                MPI_Isend(grid.data() + z_size * y_size, layer_size, MPI_DOUBLE, comm_rank - 1, 1, MPI_COMM_WORLD,&request_send_left);
            }

            if (comm_rank == RANK_ROOT) {
                MPI_Irecv(right_recv_border, layer_size, MPI_DOUBLE, comm_rank + 1, 1,MPI_COMM_WORLD, &request_recv_right);
            } else if (comm_rank == comm_size - 1) {
                MPI_Irecv(grid.data(), layer_size, MPI_DOUBLE, comm_rank - 1, 1, MPI_COMM_WORLD, &request_recv_left);
            } else {
                MPI_Irecv(right_recv_border, layer_size, MPI_DOUBLE, comm_rank + 1, 1, MPI_COMM_WORLD, &request_recv_right);
                MPI_Irecv(grid.data(), layer_size, MPI_DOUBLE, comm_rank - 1, 1, MPI_COMM_WORLD, &request_recv_left);
            }
        }

        for (int i = 2; i < x_size - 2; ++i) {
            for (int k = 0; k < z_size; ++k) {
                for (int j = 0; j < y_size; ++j) {
                    count_next_phi(i, j, k);

                }
            }
        }
        MPI_Status status;
        if (comm_size > 1) {

            if (comm_rank == RANK_ROOT) {
                MPI_Wait(&request_recv_right, &status);
            } else if (comm_rank == comm_size - 1) {
                MPI_Wait(&request_recv_left, &status);
            } else {
                MPI_Wait(&request_recv_right, &status);
                MPI_Wait(&request_recv_left, &status);
            }
            std::cout << "aftwait";
            std::cout.flush();
        }

        for (int k = 0; k < z_size; ++k) {
            for (int j = 0; j < y_size; ++j) {
                count_next_phi(1, j, k);
                count_next_phi(x_size - 2, j, k);
            }
        }

        count_diff();
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&cur_maxdif, &maxdif, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        std::cout << maxdif << " ";
    }

    void count_next_phi(int i, int j, int k) {

        int cur_pos = grid_ind_from_coords(i, j, k);

        double left_x_neighbor = grid[grid_ind_from_coords(i - 1, j, k)];
        double right_x_neighbor = grid[grid_ind_from_coords(i + 1, j, k)];

        double backw_y_neighbor = (j >= (y_size - 1)) ? phi(i - 1, j + 1, k) : grid[grid_ind_from_coords(i, j + 1, k)];
        double forw_y_neighbor = (j <= 0) ? phi(i - 1, -1, k) : grid[grid_ind_from_coords(i, j - 1, k)];

        double up_z_neighbor = (k >= (z_size - 1)) ? phi(i - 1, j, k + 1) : grid[grid_ind_from_coords(i, j, k + 1)];
        double down_z_neighbor = (k <= 0) ? phi(i - 1, j, -1) : grid[grid_ind_from_coords(i, j, k - 1)];

        double v1 = (left_x_neighbor + right_x_neighbor) / (h_x * h_x) +
                    (backw_y_neighbor + forw_y_neighbor) / (h_y * h_y) +
                    (up_z_neighbor + down_z_neighbor) / (h_z * h_z);
        double phivar = phi(i - 1, j, k);
        double next_var = coef * (v1 - 6 + a * phivar);
        next_phi[cur_pos] = next_var;
    }

    void count_diff() {
        cur_maxdif = -1;
        for (int i = y_size * z_size; i < grid_size - y_size * z_size; ++i) {
            double curdif = fabs(grid[i] - next_phi[i]);
            if (curdif > cur_maxdif) {
                cur_maxdif = curdif;
            }
            grid[i] = next_phi[i];
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
                    fprintf(file, "%.4lf ", grid[grid_ind_from_coords(k, i, j)]);
                }
                fputc('\n', file);
            }
            fputc('\n', file);
        }
        std::cout << "end write" << std::endl;
    }

    void count_delta() {
        double max_delta = -1;
        for (int k = 1; k < x_size - 1; k++) {
            for (int j = 0; j < z_size; j++) {
                for (int i = 0; i < y_size; i++) {
                    double phivar = phi(k - 1, i, j);
                    double delta = fabs(grid[grid_ind_from_coords(k, i, j)] - phivar);
                    if (delta > max_delta) {
                        max_delta = delta;
                    }
                }
            }
        }
        std::cout << max_delta;
    }

};


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << size << " " << rank << std::endl;
    int x_size;
    if (N_x % size != 0) {
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 1;
    } else {
        x_size = N_x / size;
    }
    Grid_3D grid3D(x_size, N_y, N_z, 0);

    grid3D.fill_grid(size);

    int iter = 1;
    while (grid3D.maxdif > eps) {
        //   std::cout << "Maxdif from while " <<grid3D.maxdif << " ";
        std::cout << "Iteration " << iter << '\n';
        grid3D.iteration(size);
        iter++;
    }
    //   std::cout << "Maxdif after while " <<grid3D.maxdif << " ";

    if (rank == RANK_ROOT) {
        FILE *file = fopen("result.txt", "w");
        grid3D.print_grid(file);
        fclose(file);
    }
    grid3D.count_delta();
    std::cout << std::endl;


    MPI_Finalize();
    return 0;
}

//Slice x = 0
//4.25 3.50 3.25 3.50 4.25
//3.50 2.75 2.50 2.75 3.50
//3.25 2.50 2.25 2.50 3.25
//3.50 2.75 2.50 2.75 3.50
//4.25 3.50 3.25 3.50 4.25