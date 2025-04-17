#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <vector>
#include <functional>

using namespace std;


const int step_t = 1;
const int step_x = 1;

const double c = 1;
const double rho = 1;
const double pi = 3.14159265359;


double zero_func(double x, double y) {
    return 0;
}

double test1_func(double x, double y) {
    return 1;
}

double test_seminar_left(double t, double y) {
    return sin(pi*y);
}


double test_seminar_bottom(double t, double x) {
    return sin(pi*x);
}


double initial_test_seminar(double x, double y) {
    if (abs(x) < 1e-15) {
        return test_seminar_left(0, y);
    }
    if (abs(y) < 1e-15) {
        return test_seminar_bottom(0, x);
    }
    return 0;
}

double inftyNorm(vector<double> vec1, vector<double> vec2) {
    double res = 0;
    int n = vec1.size();
    for (int i = 0; i < n; i++) {
        res = max(res, abs(vec1[i] - vec2[i]));
    }
    return res;
}


void DisplayMatrix(vector<vector<double>> Matrix) {
    for (auto& row : Matrix) {
        for (auto el : row) {
            cout << el << " ";
        }
        cout << endl;
    }
}


void DisplayMatrix2(vector<vector<double>> Matrix) {
    for (int i = Matrix.size() - 1; i >= 0; i--) {
        for (int j = 0; j < Matrix[0].size(); j++) {
            cout << Matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void DisplayVector(vector<double> vec) {
    for (auto el : vec) {
        cout << el << " ";
    }
    cout << endl;
}


vector<double> Progonka(vector<vector<double>> matrix, vector<double> right) {
    // matrix[i] - вектор содержащий 3 элта a_i b_i, c_i лишние нули не содержатся
    int n = matrix.size();
    vector<double> alpha{}, beta{}, res(n, 0);

    alpha.push_back(0);
    beta.push_back(0);
    alpha.push_back(-matrix[0][2] / matrix[0][1]);
    beta.push_back(right[0] / matrix[0][1]);
    for (int i = 2; i < n; i++) {
        double znam = matrix[i - 1][1] + matrix[i - 1][0] * alpha[i - 1];
        alpha.push_back(-matrix[i - 1][2] / znam);
        beta.push_back((right[i - 1] - matrix[i - 1][0] * beta[i - 1]) / znam);
    }

    res[n - 1] = (right[n - 1] - matrix[n - 1][0] * beta[n - 1]) / (matrix[n - 1][1] + matrix[n - 1][0] * alpha[n - 1]);
    for (int k = n - 2; k >= 0; k--) {
        res[k] = res[k + 1] * alpha[k + 1] + beta[k + 1];
    }
    return res;
}




double CalculateDer_t(vector<vector<double>> M1, vector<vector<double>> M2, double tau) {
    double res = -1;
    for (int i = 0; i < M1.size(); i++) {
        for (int j = 0; j < M2.size(); j++) {
            double cur_der = abs(M2[i][j] - M1[i][j]) / tau;
            res = max(res, cur_der);
        }
    }
    return res;
}


void write_matrix(vector<vector<double>> M, string file) {
    ofstream out(file);
    if (out.is_open()) {
        for (int i = 0; i < M.size(); i++) {
            for (int j = 0; j < M[0].size(); j++) {
                out << M[i][j] << " ";
            }
            out << endl;
        }
    }
}


vector<vector<double>> SolvePoisson(
    double T, double L1, double L2, double tau, double h1, double h2,
    double(&u_left)(double t, double y), double(&u_right)(double t, double y),
    double(&u_upper)(double t, double x), double(&u_below)(double t, double x),
    double(&f)(double x, double y), double(&initial)(double x, double y), string file,
    double eps
) {
    vector<double> gridx = { 0 }, gridy = { 0 }, gridt = { 0 };
    double t = 0, x = 0, y = 0;
    while (T - t > 5e-15) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size();
    int m = gridy.size();

    vector<vector<double>> u_next, u_prev(m, vector<double>(n, 0)), u_prev_prev;

    for (int  j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            u_prev[j][i] = initial(gridx[i], gridy[j]);
        }
    }
    //DisplayMatrix(u_prev);
    vector<vector<double>> progonka_coef;
    vector<double> d;
    for (int k = 1; k < gridt.size()-1; k++) {
        t = gridt[k];
        double t_half = t + tau / 2;
        double t_next = gridt[k+1];
        u_next = {};
        u_next.push_back(vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            u_next[0][i] = u_below(t_half, gridx[i]);
        }
        //DisplayMatrix(u_next);


        for (int j = 1; j < m-1; j++) {
            progonka_coef = {};
            d = {};
            y = gridy[j];
            progonka_coef.push_back({0,1,0});
            d.push_back(u_left(t_half, y));
            for (int i = 1; i < n-1; i++) {
                x = gridx[i];
                progonka_coef.push_back({ -1/(h1*h1), 2 / tau + 2 / (h1*h1), -1 / (h1 * h1) });
                d.push_back(2/tau * u_prev[j][i] + 1/(h2*h2)*(u_prev[j+1][i] -2 * u_prev[j][i] + u_prev[j-1][i]) + f(x, y));

            }
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_right(t_half, y));

            u_next.push_back(Progonka(progonka_coef, d));

        }
        u_next.push_back(vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            u_next[m - 1][i] = u_upper(t_half, gridx[i]);
        }
        cout << "upper" << u_upper(t_half, 0.5) << endl;
        DisplayMatrix(u_next);
        cout << endl;
        // Целый шаг
        u_prev_prev = u_prev;
        u_prev = u_next;

        u_next = {};
        for (int j = 0; j < m; j++) {
            u_next.push_back(vector<double>(n, 0));
        }
        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][0] = u_left(t_next, gridy[j]);
        }
        //DisplayMatrix(u_next);
        for (int i = 1; i < n - 1; i++) {
            progonka_coef = {};
            d = {};
            x = gridx[i];
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_below(t_next, x));
            for (int j = 1; j < m - 1; j++) {
                y = gridy[j];
                progonka_coef.push_back({ -1 / (h2 * h2), 2 / tau + 2 / (h2 * h2), -1 / (h2 * h2) });
                d.push_back(2/tau *u_prev[j][i] + 1/(h1*h1)*(u_prev[j][i+1] - 2*u_prev[j][i] + u_prev[j][i-1]) + f(x,y));
            }
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_upper(t_next, x));
            vector<double> result_progonka = Progonka(progonka_coef, d);

            //DisplayVector(result_progonka);
            for (int j = 0; j < m; j++) {
                u_next[j][i] = result_progonka[j];
            }
            //u_next.push_back();

        }

        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][n-1] = u_right(t_next, gridy[j]);
        }

        if (CalculateDer_t(u_next, u_prev_prev, tau) < eps) {
            cout << "Final t = " << t_next << endl;
            write_matrix(u_next, file);
            return u_next;
        }
        u_prev = u_next;
        
    }

}



int main()
{
    double T = 100000, L1 = 1, L2 = 1, tau = 0.1, h1 = 0.1, h2 = 0.1;
    //auto u_left = test1_func, u_right = test1_func, u_upper = test1_func, u_below = test1_func;
    vector<vector<double>> matrix = SolvePoisson(T, L1, L2, tau, h1, h2, test_seminar_left, zero_func, zero_func, test_seminar_bottom, zero_func, initial_test_seminar, "test.txt", 1e-4);
    DisplayMatrix(matrix);
}


