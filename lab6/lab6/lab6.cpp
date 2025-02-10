#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
using namespace std;
const double epsilon = 1e-8;
const int Maxiter = 1000;

double f1_system1_book(vector<double> values) {
    return 2 * values[1] + values[2] * values[2] - 1;
}


double f2_system1_book(vector<double> values) {
    return 6 * values[1] - values[2] * values[2] + 1;
}


double f1_system1_test(vector<double> values) {
    return values[2];
}


double f2_system1_test(vector<double> values) {
    return values[1];
}


void DisplayVector(vector<double> vec) {
    for (auto el : vec) {
        cout << el << " ";
    }
    cout << endl;
}


void DisplayMatrix(vector<vector<double>> Matrix) {
    for (auto& row : Matrix) {
        for (auto el : row) {
            cout << el << " ";
        }
        cout << endl;
    }
}


void triangularize(vector<vector<double>>& matrix) {
    int n = matrix.size();

    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            double temp = matrix[i][j];
            for (int k = 0; k < n + 1; ++k) {
                matrix[i][k] -= matrix[j][k] * temp / matrix[j][j];
            }
        }
    }
}


bool check_matrix(const vector<vector<double>>& triang_matrix) {
    double num = 1;
    for (int i = 0; i < triang_matrix.size(); ++i) {
        num *= triang_matrix[i][i];
    }
    return (abs(num) > 1e-15);
}


vector<double> Gauss_method(vector<vector<double>>& matrix) {
    triangularize(matrix);
    int n = matrix.size();
    vector<double> x(n);
    if (!check_matrix(matrix)) {
        cout << "det A = 0" << endl;
        return x;
    }

    double sum;
    for (int i = n - 1; i > -1; --i) {
        sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += matrix[i][j] * x[j];
        }
        x[i] = (matrix[i][n] - sum) / matrix[i][i];
    }
    return x;
}

double inftyNorm(vector<double> vec1, vector<double> vec2) {
    double res = -9999999999999999999;
    int n = vec1.size();
    for (int i = 0; i < n; i++) {
        res = max(res, abs(vec1[i] - vec2[i]));
    }
    return res;
}

vector<vector<double>> ComputeJacoby(vector<double(*)(vector<double>)> f, vector<double> Point) {
    int n = f.size();
    vector<double> NewPoint{}, str{};
    vector<vector<double>> MatrixJacoby{};
    for (int fun= 0; fun < n; fun++) {
        str = {};
        for (int var = 0; var < n; var++) {
            NewPoint = Point;
            NewPoint[var] += epsilon;
            str.push_back((f[fun](NewPoint) - f[fun](Point)) / epsilon);
        }
        MatrixJacoby.push_back(str);
    }
    return MatrixJacoby;
}


vector<double>NewtonMethod(vector<double> y0, vector<double(*)(vector<double>)> f) {

    int iteration = 0, n = f.size();
    vector<double> prevstep = y0;
    vector<double> curstep = y0, b = y0, delta{};
    vector<vector<double>> JacobyMatrix;
    while (iteration < Maxiter) {
        for (int i = 0; i < n; i++) {
            b[i] = -f[i](prevstep);
        }
        JacobyMatrix = ComputeJacoby(f, prevstep);
        // Сводим к виду чтобы скормить Гауссу
        for (int i = 0; i < n; i++) {
            JacobyMatrix[i].push_back(b[i]);
        }
        delta = Gauss_method(JacobyMatrix);
        for (int i = 0; i < n; i++) {
            curstep[i] = prevstep[i] + delta[i];
        }
        if (inftyNorm(curstep, prevstep) < 1e-6) { break; }
        prevstep = curstep;
        iteration++;
    }
    return curstep;
}

// tau - шаг, T = tmax, y0 - начанльная точка, f - вектор правых частей
vector<double> ImplicitEuler(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
    int n = f.size();
    // Делаем сетку
    vector<double> time = {};
    double t0 = 0;
    while (t0 <= T) {
        time.push_back(t0);
        t0 += tau;
    }
    //Нужно найти все значения в узлах
    vector<vector<double>> res = { y0 };
    for (int i = 1; i < time.size(); i++) {
        // Нужно найти y_i. Для этого используем метод Ньютона
        int iteration = 0;
        vector<double> prevstep = res[i-1], curstep = res[i - 1], b = res[i - 1], delta{};
        vector<vector<double>> JacobyMatrix{};

        while (iteration < Maxiter) {
            for (int j = 0; j < n; j++) {
                b[j] = -tau*f[j](time[i], prevstep)+ prevstep[j] - res[i-1][j];
            }
            vector<double> NewPoint{}, str{};
            for (int fun = 0; fun < n; fun++) {
                str = {};
                for (int var = 0; var < n; var++) {
                    NewPoint = prevstep;
                    NewPoint[var] += epsilon;
                    str.push_back((tau*f[fun](time[i],NewPoint) - f[fun](time[i],prevstep)) / epsilon);
                }
                JacobyMatrix.push_back(str);
            }
            delta = Gauss_method(JacobyMatrix);
        }
    }
}


vector<double> GenerateUniformGrid(double t0, double T, int n) {
    // n - количество отрезков, [t0, T] - отрезок 
    vector<double> res;
    double h = (T - t0) / n;
    vector<double> values{};

    for (int i = 0; i < n + 1; i++) {
        res.push_back(t0 + i * h);
    }

    return res;
}


vector<vector<double>> Euler_explicit(
    vector<double> y0,
    const vector<double>& grid,
    vector<double(*)(vector<double>)> syst) {

    int n = grid.size(); // количество узлов
    int num_vars = syst.size(); // количество переменных
    double tau = grid[1] - grid[0];
    vector<vector<double>> res; // res хранит вектор по каждой слою времени, 
    // в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i

    res.emplace_back(y0); // добавляем начальный слой времени

    // i - номер узла, j - номер переменной
    for (int i = 1; i < n; i++) {
        vector<double> y_cur(num_vars);
        vector<double> prev_vars;
        prev_vars.emplace_back(grid[i - 1]);
        prev_vars.insert(prev_vars.end(), res[i - 1].begin(), res[i - 1].end());
        for (int j = 0; j < num_vars; j++) {
            y_cur[j] = res[i - 1][j] + tau * syst[j](prev_vars);

        }
        res.emplace_back(y_cur);
    }

    return res;
}





vector<double> GenerateUniformGrid(double t0, double T, int n) {
    // n - количество отрезков, [t0, T] - отрезок 
    vector<double> res;
    double h = (T - t0) / n;
    vector<double> values{};

    for (int i = 0; i < n + 1; i++) {
        res.push_back(t0 + i * h);
    }

    return res;
}


vector<vector<double>> Euler_explicit(
    vector<double> y0,
    const vector<double>& grid,
    vector<double(*)(vector<double>)> syst) {

    int n = grid.size(); // количество узлов
    int num_vars = syst.size(); // количество переменных
    double tau = grid[1] - grid[0];
    vector<vector<double>> res; // res хранит вектор по каждой слою времени, 
    // в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i

    res.emplace_back(y0); // добавляем начальный слой времени

    // i - номер узла, j - номер переменной
    for (int i = 1; i < n; i++) {
        vector<double> y_cur(num_vars);
        vector<double> prev_vars;
        prev_vars.emplace_back(grid[i - 1]);
        prev_vars.insert(prev_vars.end(), res[i - 1].begin(), res[i - 1].end());
        for (int j = 0; j < num_vars; j++) {
            y_cur[j] = res[i - 1][j] + tau * syst[j](prev_vars);

        }
        res.emplace_back(y_cur);
    }

    return res;
}





void PrintGridFunc(const vector<vector<double>>& vec) {
    int num_vars = vec[0].size();
    int n = vec.size();
    cout << endl;
    for (int j = 0; j < num_vars; j++) {
        cout << "x" << j + 1 << " ";
    }
    cout << endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < num_vars; j++) {
            cout << vec[i][j] << " ";
        }
        cout << endl;
    }
}

int main()
{
    
}
