#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <iomanip>


using namespace std;
const int M = 160;//40,80
const int N = 180;//40,90
const double xmin = -3.0, xmax = 3.0;
const double ymin = 0.0, ymax = 4.0;
const double hx = (xmax - xmin) / M;
const double hy = (ymax - ymin) / N;
const double tolerance = 2e-4; 
const double epsilon = 1e-2;
double norm = 1.0;
int nx =M + 1 ;
int ny =N + 1;
int iters=1;
int test = 1;

bool isInside(double x, double y) {
    return  (abs(4 * x) + 3 * y - 12 <= 0)&& y>=0;
}

double ly(double x, double y) {
    double x_r=x - 0.5 * hx;
    double y_up=y + 0.5 * hy,y_dw= y - 0.5 * hy;

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

double ls(double x, double y) {
    double x_l = x - 0.5 * hx;
    double y_dw = y - 0.5 * hy;
    double x_r = x + 0.5 * hx;
    double y_up = y + 0.5 * hy;

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


double lx(double x, double y) {
    double y_dw= y - 0.5 * hy;
    double x_l=x - 0.5 * hx,x_r= x + 0.5 * hx;

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


vector<vector<double>> speed(vector<vector<double>>& u, vector<vector<double>>& A, vector<vector<double>>& B, vector<vector<double>>& F) {
    while (norm > tolerance) {
        vector<vector<double>> r(M + 1, vector<double>(N + 1, 0.0));
        vector<vector<double>> Ar(M + 1, vector<double>(N + 1, 0.0));
        norm = 0.0;
        double alpha;
        double R = 0.0, AR = 0.0;

#pragma omp parallel
        {
            double local_R = 0.0;
            double local_AR = 0.0;

#pragma omp for
            for (int i = 1; i < M; i++) {
                for (int j = 1; j < N; j++) {
                    double x = xmin + i * hx;
                    double y = ymin + j * hy;
                    if (isInside(x, y)) {
                        r[i][j] = -F[i][j] - 1.0 / hx * (A[i + 1][j] * (u[i + 1][j] - u[i][j]) / hx - A[i][j] * (u[i][j] - u[i - 1][j]) / hx)
                            - 1.0 / hy * (B[i][j + 1] * (u[i][j + 1] - u[i][j]) / hy - B[i][j] * (u[i][j] - u[i][j - 1]) / hy);
                    }
                    else {
                        r[i][j] = 0.0;
                    }
                }
            }

#pragma omp for
            for (int i = 1; i < M; i++) {
                for (int j = 1; j < N; j++) {
                    double x = xmin + i * hx;
                    double y = ymin + j * hy;
                    if (isInside(x, y)) {
                        Ar[i][j] = -1.0 / hx * (A[i + 1][j] * (r[i + 1][j] - r[i][j]) / hx - A[i][j] * (r[i][j] - r[i - 1][j]) / hx)
                            - 1.0 / hy * (B[i][j + 1] * (r[i][j + 1] - r[i][j]) / hy - B[i][j] * (r[i][j] - r[i][j - 1]) / hy);
                    }
                    else {
                        Ar[i][j] = 0.0;
                    }
                }
            }

#pragma omp for
            for (int i = 1; i < M; i++) {
                for (int j = 1; j < N; j++) {
                    local_R += r[i][j] * r[i][j];
                    local_AR += Ar[i][j] * r[i][j];
                }
            }

#pragma omp critical
            {
                R += local_R;
                AR += local_AR;
            }
        }

        alpha = R / AR;

        for (int i = 1; i < M; i++) {
            for (int j = 1; j < N; j++) {
                u[i][j] = u[i][j] - alpha * r[i][j];
                norm += r[i][j] * r[i][j];
            }
        }

        norm = sqrt(norm) * alpha;
        if (norm < tolerance) {
            break;
        }
        iters++;
    }
    return u;
}
void saveToCSV(const vector<vector<double>>& u, const string& filename) {
    ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < M+1; ++i) {
            for (int j = 0; j < N+1; ++j) {
                file << u[i][j];
                if (j < N) {
                    file << ","; 
                }
            }
            file << endl; 
        }
        file.close();
        cout << "Matrix saved to " << filename << endl;
    }
    else {
        cerr << "Failed to open file: " << filename << endl;
    }
}


int main() {
    vector<vector<double>> u(nx, vector<double>(ny, 0.0)); 
    vector<vector<double>> A(nx, vector<double>(ny, 0.0)); 
    vector<vector<double>> B(nx, vector<double>(ny, 0.0)); 
    vector<vector<double>> F(nx, vector<double>(ny, 0.0)); 
    vector<vector<double>> w(nx, vector<double>(ny, 0.0));

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = xmin + i * hx;
            double y = ymin + j * hy;
            A[i][j] = ly(x, y) / hy + (1 - ly(x, y) / hy) / epsilon;
            B[i][j] = lx(x, y) / hx + (1 - lx(x, y) / hx) / epsilon;
            F[i][j] = ls(x,y)/(hx*hy);    

        }
    }
   
    auto start = chrono::high_resolution_clock::now();
    u=speed(u,A,B,F);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "M:  " << M << endl;
    cout << "N:  " << N << endl;
    cout << "iter:  " << iters << endl;
    cout << "Runtime:  " << duration.count() << "  ms" << endl;
    saveToCSV(u, "resultfinal.csv");
    return 0;
}