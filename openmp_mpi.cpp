#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include <fstream>
#include <iomanip>


using namespace std;
const int M = 40;//80,160
const int N = 40;//90,180
const double xmin = -3.0, xmax = 3.0;
const double ymin = 0.0, ymax = 4.0;
const double hx = (xmax - xmin) /M;
const double hy = (ymax - ymin) /N;
const double tolerance = 2e-4; 
const double epsilon = 1e-2;


bool isInside(double x, double y) {
    return  (abs(4 * x) + 3 * y - 12 <= 0)&& y>=0;
}

double ly(int i, int j) {
    double x_r=xmin +(i- 0.5) * hx;
    double y_up=ymin+(j + 0.5) * hy,y_dw= ymin +(j- 0.5) * hy;

    bool topInside = isInside(x_r, y_up);
    bool bottomInside =isInside(x_r, y_dw);
    
    if (topInside && bottomInside) return hy;
    if (!topInside && !bottomInside) return 0;
    
    if (topInside && !bottomInside){
        return abs(y_up) ;
    }
     if (!topInside && bottomInside) {
        return x_r>=0 ? abs(((-4)*x_r+12)/3 -y_dw): abs((4*x_r+12)/3 -y_dw);
    }

}

double ls(int i, int j) {
    double x_l = xmin +(i- 0.5) * hx;
    double y_dw = ymin+ (j- 0.5) * hy;
    double x_r = xmin +(i + 0.5)  * hx;
    double y_up = ymin + (j+0.5) * hy;

    bool l_topInside = isInside(x_l, y_up);
    bool l_bottomInside = isInside(x_l, y_dw);
    bool r_topInside = isInside(x_r, y_up);
    bool r_bottomInside = isInside(x_r, y_dw);

    double area_full = hx * hy;
    double pow_expr1 = pow((3.0 - 3.0 / 4.0 * y_dw - x_l), 2) * 4.0 / 3.0;
    double pow_expr2 = pow((3.0 - 3.0 / 4.0 * y_dw + x_r), 2) * 4.0 / 3.0;
    double pow_expr3 = pow((3.0 - 3.0 / 4.0 * y_up - x_r), 2) * 4.0 / 3.0;
    double pow_expr4 = pow((3.0 - 3.0 / 4.0 * y_up + x_l), 2) * 4.0 / 3.0;
    double mid_y = (y_dw + y_up) * 3.0 / 4.0;

    if (l_topInside && l_bottomInside && r_topInside && r_bottomInside) {
        return area_full;
    } 
    else if (l_topInside && !l_bottomInside && r_topInside && !r_bottomInside) {
        return area_full / 2.0;
    } 
    else if (!l_topInside && l_bottomInside && !r_topInside && r_bottomInside) {
        return area_full / 4.0;
    } 
    else if (!l_topInside && l_bottomInside && !r_topInside && !r_bottomInside) {
        return pow_expr1;
    } 
    else if (!l_topInside && !l_bottomInside && !r_topInside && r_bottomInside) {
        return pow_expr2;
    } 
    else if (l_topInside && l_bottomInside && !r_topInside && r_bottomInside) {
        return area_full - pow_expr3;
    } 
    else if (l_topInside && !l_bottomInside && !r_topInside && r_bottomInside) {
        return area_full - pow_expr4;
    } 
    else if (l_topInside && l_bottomInside && !r_topInside && !r_bottomInside) {
        return (mid_y + x_r - x_l) * hy / 2.0;
    } 
    else if (!l_topInside && !l_bottomInside && r_topInside && r_bottomInside) {
        return (mid_y - x_r + x_l) * hy / 2.0;
    } 
    return 0;
}


double lx(int i, int j) {
    double y_dw= ymin +(j- 0.5) * hy;
    double x_l=xmin + (i- 0.5) * hx,x_r= xmin+ (i+ 0.5) * hx;

    bool leftInside = isInside(x_l, y_dw);
    bool rightInside =isInside(x_r, y_dw);
    
    if (leftInside && rightInside) return hx;
    if (!leftInside && !rightInside) return 0;
    
    if (leftInside && !rightInside) {
        return  abs((12-3*y_dw)/4-x_l) ;
    }
     if (!leftInside && rightInside) {
        return  abs(x_r - (3*y_dw-12)/4); 
    }

}

// 动态分割矩形域
void decompose_domain(int rank, int size, int M, int N, int &start_row, int &end_row, int &start_col, int &end_col) {
    int Px = 1, Py = size;
    for (int i = 1; i <= size; ++i) {
        if (size % i == 0) {
            int potential_Px = i;
            int potential_Py = size / i;
            double ratio = static_cast<double>(potential_Px) / potential_Py;
            if (ratio >= 0.5 && ratio <= 2.0) {
                Px = potential_Px;
                Py = potential_Py;
            }
        }
    }

    int x = rank % Px;
    int y = rank / Px;

    int rows_per_domain = M / Py;
    int row_remainder = M % Py;
    start_row = y * rows_per_domain + min(y, row_remainder);
    end_row = start_row + rows_per_domain + (y < row_remainder ? 1 : 0) - 1;

    int cols_per_domain = N / Px;
    int col_remainder = N % Px;
    start_col = x * cols_per_domain + min(x, col_remainder);
    end_col = start_col + cols_per_domain + (x < col_remainder ? 1 : 0) - 1;
}

// 并行最速下降法
void computeUij(vector<vector<double>> &u, const vector<vector<double>> &A, const vector<vector<double>> &B, const vector<vector<double>> &F,
                int M, int N, double hx, double hy, double tolerance,
                int rank, int neighbor_top, int neighbor_bottom, int neighbor_left, int neighbor_right,int &iterations) {
    double norm = 1.0;
    iterations = 0;

    vector<double> send_top(N, 0.0), send_bottom(N, 0.0), send_left(M, 0.0), send_right(M, 0.0);
    vector<double> recv_top(N, 0.0), recv_bottom(N, 0.0), recv_left(M, 0.0), recv_right(M, 0.0);

    while (norm > tolerance) {
        vector<vector<double>>  r(M, vector<double>(N, 0.0));
        vector<vector<double>>  Ar (M, vector<double>(N, 0.0));
        double local_R = 0.0, local_AR = 0.0, R = 0.0, AR = 0.0;
        #pragma omp parallel for
        for (int j = 1; j < N - 1; ++j) {
            send_top[j] = u[1][j];
            send_bottom[j] = u[M - 2][j];
        }
        #pragma omp parallel for
        for (int i = 1; i < M - 1; ++i) {
            send_left[i] = u[i][1];
            send_right[i] = u[i][N - 2];
        }

        MPI_Sendrecv(send_top.data(), N, MPI_DOUBLE, neighbor_top, 0,
                     recv_bottom.data(), N, MPI_DOUBLE, neighbor_bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_bottom.data(), N, MPI_DOUBLE, neighbor_bottom, 1,
                     recv_top.data(), N, MPI_DOUBLE, neighbor_top, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_left.data(), M, MPI_DOUBLE, neighbor_left, 2,
                     recv_right.data(), M, MPI_DOUBLE, neighbor_right, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_right.data(), M, MPI_DOUBLE, neighbor_right, 3,
                     recv_left.data(), M, MPI_DOUBLE, neighbor_left, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (neighbor_top >= 0)
            #pragma omp parallel for
            for (int j = 1; j < N - 1; ++j) u[0][j] = recv_top[j];
        if (neighbor_bottom >= 0)
            #pragma omp parallel for
            for (int j = 1; j < N - 1; ++j) u[M - 1][j] = recv_bottom[j];
        if (neighbor_left >= 0)
            #pragma omp parallel for
            for (int i = 1; i < M - 1; ++i) u[i][0] = recv_left[i];
        if (neighbor_right >= 0)
            #pragma omp parallel for
            for (int i = 1; i < M - 1; ++i) u[i][N - 1] = recv_right[i];
        #pragma omp parallel for reduction(+:local_R, local_AR)
        for (int i = 1; i < M - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                r[i][j] = -F[i][j]
                          - 1.0 / hx * (A[i + 1][j] * (u[i + 1][j] - u[i][j]) / hx - A[i][j] * (u[i][j] - u[i - 1][j]) / hx)
                          - 1.0 / hy * (B[i][j + 1] * (u[i][j + 1] - u[i][j]) / hy - B[i][j] * (u[i][j] - u[i][j - 1]) / hy);

                Ar[i][j] = -1.0 / hx * (A[i + 1][j] * (r[i + 1][j] - r[i][j]) / hx - A[i][j] * (r[i][j] - r[i - 1][j]) / hx)
                           - 1.0 / hy * (B[i][j + 1] * (r[i][j + 1] - r[i][j]) / hy - B[i][j] * (r[i][j] - r[i][j - 1]) / hy);

                local_R += r[i][j] * r[i][j];
                local_AR += Ar[i][j] * r[i][j];
            }
        }

        MPI_Allreduce(&local_R, &R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_AR, &AR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double alpha = R / AR;
        norm = 0.0;
        #pragma omp parallel for reduction(+:norm)
        for (int i = 1; i < M - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                u[i][j] -= alpha * r[i][j];
                norm += r[i][j] * r[i][j];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        norm = sqrt(norm) * alpha;

        iterations++;
    }
}


int main(int argc, char *argv[]) {
    //MPI初始化
    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//总进程数
    MPI_Comm_size(MPI_COMM_WORLD, &size);//进程编号

    // 子域划分
    int start_row, end_row, start_col, end_col;
    decompose_domain(rank, size, M, N, start_row, end_row, start_col, end_col);

    // 子域大小
    int local_M = end_row - start_row + 3; // 包括幽灵节点
    int local_N = end_col - start_col + 3;

    // 初始化子域矩阵
    vector<vector<double>>  u(local_M, vector<double>(local_N, 0.0));
    vector<vector<double>>  A(local_M, vector<double>(local_N, 0.0));
    vector<vector<double>>  B(local_M, vector<double>(local_N, 0.0));
    vector<vector<double>>  F(local_M, vector<double>(local_N, 0.0));

    // 计算 A, B, F 矩阵
    #pragma omp parallel for
    for (int i = start_row; i <= end_row; ++i) {
        for (int j = start_col; j <= end_col; ++j) {
            
            A[i - start_row + 1][j - start_col + 1] = ly(i, j) / hy + (1 - ly(i, j) / hy) / epsilon;;
            B[i - start_row + 1][j - start_col + 1] = lx(i, j) / hx + (1 - lx(i, j) / hx) / epsilon;;
            F[i - start_row + 1][j - start_col + 1] = ls(i, j) / (hx * hy);
        }
    }

    // 邻居进程编号
    int neighbor_top = (start_row > 0) ? rank - size / 2 : MPI_PROC_NULL;
    int neighbor_bottom = (end_row < M - 1) ? rank + size / 2 : MPI_PROC_NULL;
    int neighbor_left = (start_col > 0) ? rank - 1 : MPI_PROC_NULL;
    int neighbor_right = (end_col < N - 1) ? rank + 1 : MPI_PROC_NULL;

    // 计算 u 矩阵
    int iterations = 0;
    computeUij(u, A, B, F, local_M, local_N, hx, hy, tolerance, rank,
               neighbor_top, neighbor_bottom, neighbor_left, neighbor_right,iterations);

    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    int max_iterations;
    MPI_Reduce(&iterations, &max_iterations, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // 主进程保存全局结果并输出迭代次数和运行时间
    if (rank == 0) {
       
        cout << "Total iterations: " << max_iterations << endl;
        cout << "Total runtime: " << max_time << " seconds" << endl;
    }


    MPI_Finalize();
    return 0;
}