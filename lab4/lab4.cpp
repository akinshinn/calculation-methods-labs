
#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include "Polynom.h"
# define pi 3.14159265358979323846
const int koltestpoints = 1000;
using namespace std;
void DisplayMatrix(vector<vector<double>> Matrix) {
    for (auto& row : Matrix) {
        for (auto el : row) {
            cout << el << " ";
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
    vector<double> alpha{}, beta{}, res(n,0);

    alpha.push_back(0);
    beta.push_back(0);
    alpha.push_back(-matrix[0][2] / matrix[0][1]);
    beta.push_back(right[0] / matrix[0][1]);
    for (int i = 2; i < n ; i++) {
        double znam = matrix[i - 1][1] + matrix[i - 1][0] * alpha[i - 1];
        alpha.push_back(-matrix[i - 1][2] / znam);
        beta.push_back((right[i - 1] - matrix[i - 1][0] * beta[i - 1]) / znam);
    }
    
    res[n-1] = (right[n - 1] - matrix[n - 1][0] * beta[n - 1]) / (matrix[n - 1][1] + matrix[n - 1][0] * alpha[n - 1]);
    for (int k = n - 2; k >=0; k--) {
        res[k] = res[k + 1] * alpha[k + 1] + beta[k + 1];
    }
    return res;
}




double func0(double x) {
    return x;
}

double func1(double x) {
    return x * x;
}

double func2(double x) {
    return 1 / (1 + x * x);
}

double func3(double x) {
    return 1 / atan(1 + 10 * x * x);
}

double func4(double x) {
    return pow((4 * x * x * x + 2 * x * x - 4 * x + 2), sqrt(2)) + asin(1 / (5 + x - x * x)) - 5;
}
double func5(double x) {
    return exp(x);
}

double func6(double x) {
    return sin(pi * x);
}
double runge(double x) {
    return 1 / (1 + 25*x * x);
}

double func_test(double x) {
    double R1 = cos(pow(x,5) - x + 3 + cbrt(2));
    double R2 = atan((pow(x,3) - 5 * sqrt(2) * x - 4) / (sqrt(6) * x + sqrt(3)));

    return R1 + R2 + 1.8;
}

double const_f(double x) {
    return 1;
}


pair<vector<double>, vector<double>> GenerateUniformGrid(double a, double b, int n, double(& func)(double x)) {
    // n - количество узлов, [a, b] - отрезок интерполирования
    vector<double> xValues{}, yValues{};
    double h = (b - a) / n;
    for (int i = 0; i < n + 1; i++) {
        xValues.push_back(a + i * h);
        yValues.push_back(func(a + i * h));
    }
    return { xValues, yValues };
}



pair<vector<double>, vector<double>> GenerateChebyshevGrid(double a, double b, int n, double(&func)(double x)) {
    vector<double> xValues{}, yValues{};
    for (int i = 0; i < n + 1; i++) {
        double xcurr = (a + b) / 2 + (b - a) / 2 * cos(((2 * i + 1) * pi) / (2 * (n + 1)));
        xValues.push_back(xcurr);
        yValues.push_back(func(xcurr));
    }
    return { xValues, yValues };
}

void PrintGrid(vector<double>xValues, vector<double> yValues) {
    for (int i = 0; i < xValues.size(); i++) {
        cout << xValues[i] << " " << yValues[i] << endl;
    }
}


void PrintGrid(pair<vector<double>, vector<double>> grid) {
    vector<double>& xValues = grid.first;
    vector<double>& yValues = grid.second;
    for (int i = 0; i < xValues.size(); i++) {
        cout << xValues[i] << " " << yValues[i] << endl;
    }
}


vector<vector<double>> SplainInterpolationUniformGrid(double a, double b, int n, double(&func)(double x)) {
    pair<vector<double>, vector<double>> grid = GenerateUniformGrid(a, b, n, func);
    vector<double> xValues = grid.first, yValues = grid.second;
    //cout << "Uniform grid" << endl;
    //PrintGrid(grid.first, grid.second);
    vector<double> A(n+1, 0), B(n+1, 0), C(n+1, 0), D(n+1, 0), H(n+1,0), G(n+1,0);
    
    for (int i = 1; i < n+1; i++) {
        A[i] = yValues[i - 1];
        H[i] = xValues[i] - xValues[i - 1];
        G[i] = (yValues[i] - yValues[i - 1]) / H[i];

    }
    vector<vector<double>> SystemMatrix{};
    vector <double> right{}, sol{};
    SystemMatrix.push_back({ 0,2 * (H[1] + H[2]), H[2] });
    right.push_back(3 * (G[2] - G[1]));
    for (int i = 3; i < n; i++) {
        SystemMatrix.push_back({ H[i - 1],2 * (H[i - 1] + H[i]),H[i] });
        right.push_back(3*(G[i]-G[i-1]));
    }
    SystemMatrix.push_back({ H[n - 1],2 * (H[n - 1] + H[n]),H[n] });
    right.push_back(3 * (G[n] - G[n - 1]));
    sol = Progonka(SystemMatrix, right);

    for (int i = 2; i < n+1; i++) {
        C[i] = sol[i - 2];
    }

    for (int i = 1; i < n; i++) {
        B[i] = G[i] - (C[i + 1] + 2 * C[i]) * H[i] / 3;
        D[i] = (C[i + 1] - C[i]) / (3 * H[i]);
    }
    B[n] = G[n] - 2 * C[n] * H[n] / 3;
    D[n] = -C[n] / (3 * H[n]);

    return { A,B,C,D };
}


double split_diff(const vector<double>& x, const vector<double>& y) {
    double res = 0;
    int k = x.size();
    for (int j = 0; j < k; j++) {
        double prod = 1;
        for (int i = 0; i < j; i++) {
            prod *= x[j] - x[i];
        }
        for (int i = j+1; i < k; i++) {
            prod *= x[j] - x[i];
        }
        res += y[j] / prod;
    }
    return res;
}
double split_diff(const vector<double>& x, const vector<double>& y, int k) {
    double res = 0;
    for (int j = 0; j < k+1; j++) {
        double prod = 1;
        for (int i = 0; i < j; i++) {
            prod *= x[j] - x[i];
        }
        for (int i = j + 1; i < k+1; i++) {
            prod *= x[j] - x[i];
        }
        res +=  y[j]/ prod;
    }
    return res;
}


Polynom LagrangeInterpolation(const pair<vector<double>, vector<double>>& grid) {
    int n = grid.first.size()-1;
    const vector<double>& x = grid.first;
    const vector<double>& y = grid.second;
    Polynom p(n);
    p.coefs[n] = y[0];
    Polynom subp(vector<double> {1});
    Polynom coef;
    for (int i = 1; i <= n; i++) {
        Polynom curP({ 1,-x[i-1] });
        subp *= curP;
        coef = subp;
        subp *= split_diff(x, y, i);
        p += subp;

        subp = coef;
    }
    return p;
}


bool check_interpolate(double a, double b, double(&func)(double x), double eps, Polynom p) {
    pair<vector<double>, vector<double>> subgrid = GenerateUniformGrid(a, b, 1000, func);
    double maxx = -1;
    double value;
    for (int i = 0; i < subgrid.first.size(); i++) {
        value = abs(p.getValue(subgrid.first[i]) - subgrid.second[i]);
        if (value > eps) { 
            cout << value << endl;
            cout << "x = " << subgrid.first[i] << endl;
            cout << "y = " << subgrid.second[i] << endl;
            return true; 
        };
    }
    return false;
}


double error_norm(const Polynom& Ln, double(&func)(double x), double a, double b) {
    int n = 1000;
    double maxx = -1;
    pair<vector<double>, vector<double>> subgrid = GenerateUniformGrid(a, b, 1000, func);
    for (int i = 0; i < n; i++) {
        maxx = max(maxx, abs(subgrid.second[i] - Ln.getValue(subgrid.first[i])));
    }
    return maxx;
}


Polynom find_best_chebyshev(double a, double b, double(&func)(double x), double eps) {
    Polynom p;
    int n = 1;
    do {
        n++;
        cout << n << endl;
        pair<vector<double>, vector<double>> grid = GenerateChebyshevGrid(a, b, n, func);
        p = LagrangeInterpolation(grid);
    } while (check_interpolate(a, b, func, eps, p));
    return p;

}

double calculateNorm(double a, double b, double h, double(&func)(double x)) {
    int n = (b - a) / h;
    vector<double> A, B, C, D;
    vector<vector<double>> res = SplainInterpolationUniformGrid(a, b, n, func);
    A = res[0];
    B = res[1];
    C = res[2];
    D = res[3];
    double curr = a;
    double norm = 0;
    vector<double> xValues{}, yValues{};

    for (int i = 0; i < n + 1; i++) {
        xValues.push_back(a + i * h);
    }

    for (int i = 0; i < n; i++) {
        curr = a;
        for (int j = 0; j < 50; j++) {
            norm = max(norm, abs(func(curr) - (A[i] + B[i] * (curr - xValues[i]) + C[i] * pow(curr - xValues[i], 2) + D[i] * pow(curr - xValues[i], 3))));
            curr += h / 50;
        }
    }
    return norm;
}
int main()
{
    //double a1=-1, a2=-1, a3=-3, a4=-1, b1=1, b2=1, b3=3, b4=1;
    //int n1=5, n2=5, n3=5, n4=5;

    //vector<vector<double>> ex1 = SplainInterpolationUniformGrid(a1, b1, n1, func2);
    //cout << "Splaine coefs with Uniform Grid" << endl;
    //for (int i = 1; i <= n1; i++) {
    //    cout << ex1[0][i] << " " << ex1[1][i] << " " << ex1[2][i] << " " << ex1[3][i] << endl;
    //}
    //Polynom a((vector<double>{1}));
    //cout << a.n << endl;
    //a *= a;
    //a.print();

    //a.print();
    //a.write_file("test.txt");
    //cout << a.getValue(2);
   
<<<<<<< HEAD
    int n = 16;
    pair<vector<double>, vector<double>> grid = GenerateUniformGrid(-1,1,n,func4);
    pair<vector<double>, vector<double>> chebyshev = GenerateChebyshevGrid(-1, 1, n, func4);
    PrintGrid(grid);    
    //cout << split_diff(grid.first, grid.second,4) << endl;
    Polynom p = LagrangeInterpolation(grid);
    p.write_file("test16_u.txt");
    p.print();
    Polynom p2 = LagrangeInterpolation(chebyshev);
    p2.write_file("test16_c.txt");
    cout << "Uniform grid error norm = " << error_norm(p, func4, -1, 1) << endl;
=======
    //int n = 128;
    //pair<vector<double>, vector<double>> grid = GenerateUniformGrid(-1,1,n,const_f);
    //pair<vector<double>, vector<double>> chebyshev = GenerateChebyshevGrid(-1, 1, n, const_f);
    //PrintGrid(grid);    
    //cout << split_diff(grid.first, grid.second,4) << endl;
    //Polynom p = LagrangeInterpolation(grid);
    //p.print();
    //p.write_file("constF128_u.txt");
    //Polynom p2 = LagrangeInterpolation(chebyshev);
    //p2.print();
    //p2.write_file("constF128_c.txt");
    //cout << "Uniform grid error norm = " << error_norm(p, runge, -1, 1) << endl;
>>>>>>> 7cdee46c01cf605ac895434ec51f6a4b4065ba36
    //cout << "Chebyshev grid error norm = " << error_norm(p2, runge, -1, 1) << endl;
    ////cout << p.getValue(1.5);
    //p.write_file("exp.txt");
    
    //Polynom p = find_best_chebyshev(0, 10, exp, 0.0001);
    //p.write_file("atan.txt");
    //p.print();

    double q = 0.5, h = 0.5, a = -1, b = 1;
    for (int i = 0; i < 5; i++) {
        cout << "err" << calculateNorm(a, b, pow(q, i) * h, func6) << endl;
    }
}
