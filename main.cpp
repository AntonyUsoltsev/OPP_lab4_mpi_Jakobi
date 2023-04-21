#include <iostream>
#include <vector>
#include <cmath>

constexpr double a = 1e5;
constexpr double eps = 1e-8;
constexpr double Dx = 2, Dy = 2, Dz = 2;
constexpr double N_x = 5, N_y = 5, N_z = 5;
constexpr double x_0 = -1, y_0 = -1, z_0 = -1;

double h_value(int D, int N) {
    return (double) D / (N - 1);
}

const double h_x = Dx / (N_x - 1);
const double h_y = Dy / (N_y - 1);
const double h_z = Dz / (N_z - 1);
//const double coef = 1 / (2 / (h_x * h_x) + 2 / (h_y * h_y) + 2 / (h_z * h_z) + a);
const double coef = 1 / (2 / (h_x * h_x) + a);

class Grid_1D {
public:
    std::vector<double> grid;
    std::vector<double> next_phi;
    bool stop_flag = true;

    Grid_1D() {
        grid.resize(N_x);
        next_phi.resize(N_x);
    }

    void fill_grid() {
        for (int i = 0; i < N_x; ++i) {
            grid[i] = 0;
        }
//        grid[0] = 1;
//        grid[N_x-1] = 1;
    }

    void count_next_phi(int i) {

//        double var1 = (j * N_x + i + N_x) > (N_x * N_y - 1) ? 1 : cur_phi[j * N_x + i + N_x];
//        double var2 = (j * N_x + i - N_x) < 0 ? 1 : cur_phi[j * N_x + i - N_x];
        double left_neighbor = (i - 1) < 0 ? 1 : grid[i - 1];
        double right_neighbor = (i + 1) >= N_x ? 1 : grid[i + 1];
        double coefc = coef;
        double v1 = (right_neighbor + left_neighbor) / (h_x * h_x);
        double ohivar = phi(i, 0, 0);
        double next_var = coefc * (v1 - 6 + a * ohivar);
        next_phi[i] = next_var;

    }


    double phi(int i, int j, int k) {
        return (x_0 + i * h_x) * (x_0 + i * h_x);
    }

    void iteration() {
        for (int i = 0; i < N_x; ++i) {
            count_next_phi(i);
        }
        double maxdif = 0;
        for (int i = 0; i < N_x; ++i) {
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
        for (int i = 0; i < N_x; ++i) {
            //std::cout << grid[i] << " ";
            printf("%.9lf ", grid[i]);
        }
        std::cout << std::endl;
    }

};


int main() {
    Grid_1D grid1D;
    grid1D.fill_grid();
    grid1D.print_grid();
    // std::cout << h_x << " " << N_x << " " << Dx << " " << coef<<" "<<(2 / (h_x * h_x) + a) ;
//    for (int i = 0; i < 10; i++) {
//        std::cout << "Iteration " << i+1 << '\n';
//        grid1D.iteration();
//        grid1D.print_grid();
//    }
    std::cout<<"\n\n";
    Grid_1D grid1D2;
    grid1D2.fill_grid();
int iter = 1;
    while(grid1D2.stop_flag){
        std::cout << "Iteration " << iter << '\n';
        grid1D2.iteration();
        grid1D2.print_grid();
        iter++;
    }


    return 0;
}
