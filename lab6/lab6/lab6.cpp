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
#include <functional>


const double eps = 1e-6;
using namespace std;
const double epsilon = 1e-8;
const int Maxiter = 1000;
const double k = 10;
const double m = 1;

double f2_system0_test(double t, vector<double> values) {
    return -k / m * values[0];
}

double f1_system0_test(double t, vector<double> values) {
    return values[1];
}


double f1_system1_book(double t, vector<double> values) {
    return 2 * values[0] + values[1] * values[1] - 1;
}


double f2_system1_book(double t, vector<double> values) {
    return 6 * values[0] - values[1] * values[1] + 1;
}


double f1_system1_test(double t, vector<double> values) {
    return values[1];
}


double f2_system1_test(double t, vector<double> values) {
    return values[0];
}

double f1_system2_test(double t, vector<double> values) {
    return 4 * values[0] - values[1]+exp(3*t)*(t+sin(t));
}


double f2_system2_test(double t, vector<double> values) {
    return values[0] + 2*values[1]+t*exp(3*t)*cos(t);
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
    double res = 0;
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


vector<vector<double>> ComputeJacoby(vector<function<double(vector<double>)>> f, vector<double> Point) {
    int n = f.size(); // количество функций
    vector<double> NewPoint{}, str{};
    vector<vector<double>> MatrixJacoby{};
    for (int fun = 0; fun < n; fun++) {
        str = {};
        for (int var = 1    ; var < n+1; var++) {
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



// далее везде считаем, что y0 = {y1, y2, ... yn}
vector<double>NewtonMethod(double t0, vector<double> y0, vector<function<double(vector<double>)>> f) {
    int iteration = 0, n = f.size();
    vector<double> prevstep = y0;
    vector<double> prevPoint;
    prevPoint.emplace_back(t0);
    prevPoint.insert(prevPoint.end(), prevstep.begin(), prevstep.end());
    vector<double> curstep(n), b(n), delta{};
    vector<vector<double>> JacobyMatrix;
    while (iteration < Maxiter) {
        for (int i = 0; i < n; i++) {
            b[i] = -f[i](prevPoint); // правая часть системы
        }


        JacobyMatrix = ComputeJacoby(f, prevPoint);
        // Сводим к виду чтобы скормить Гауссу
        for (int i = 0; i < n; i++) {
            JacobyMatrix[i].push_back(b[i]);
        }

        delta = Gauss_method(JacobyMatrix);


        for (int i = 0; i < n; i++) {
            curstep[i] = prevstep[i] + delta[i];
            prevPoint[i + 1] = curstep[i]; // i+1, т.к. prevPoint[0] = t0
        }

        if (inftyNorm(curstep, prevstep) < 1e-6) { break; }
        prevstep = curstep;
        iteration++;
    }
    return curstep;
}


// tau - шаг, T = tmax, y0 - начанльная точка, f - вектор правых частей
vector<vector<double>> ImplicitEuler(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
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
        vector<double> prevstep = res[i - 1], curstep = res[i - 1], b = res[i - 1], delta{};

        while (iteration < Maxiter) {
            for (int j = 0; j < n; j++) {
                b[j] = -tau * f[j](time[i], prevstep) + prevstep[j] - res[i - 1][j];
            }
            vector<vector<double>> JacobyMatrix{};
            vector<double> NewPoint = prevstep, str{};
            for (int fun = 0; fun < n; fun++) {
                str = {};
                for (int var = 0; var < n; var++) {
                    NewPoint[var] += epsilon;
                    if (fun == var) {
                        str.push_back(tau * (f[fun](time[i], NewPoint) - f[fun](time[i], prevstep)) / epsilon - 1);
                    }
                    else {
                        str.push_back(tau * (f[fun](time[i], NewPoint) - f[fun](time[i], prevstep)) / epsilon);
                    }
                    NewPoint[var] -= epsilon;
                }
                JacobyMatrix.push_back(str);
            }
            for (int j = 0; j < n; j++) {
                JacobyMatrix[j].push_back(b[j]);
            }
            delta = Gauss_method(JacobyMatrix);
            for (int j = 0; j < n; j++) {
                curstep[j] = prevstep[j] + delta[j];
            }
            if (inftyNorm(curstep, prevstep) < 1e-6) { break; }
            prevstep = curstep;
            iteration++;
        }
        res.push_back(curstep);
    }
    return res;
}

vector<vector<double>> AdamsBashford(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
    // Сетка
    int n = f.size();
    vector<double> time = {};
    double t0 = 0;
    while (t0 <= T) {
        time.push_back(t0);
        t0 += tau;
    }


    // Найдем y1, y2, y3 с помощью метода Эйлера
    vector<vector<double>> res{ y0 };

    for (int i = 1; i < 4; i++) {
        vector<double> tmp{};

        for (int j = 0; j < n; j++) {
            tmp.push_back(res[i - 1][j] + tau * f[j](time[i - 1], res[i - 1]));
        }
        res.push_back(tmp);
    }

    // Продолжим с помощью метода Адамса Бэшфорда 4-го порядка
    for (int i = 4; i < time.size(); i++) {
        vector<double> tmp{};
        for (int j = 0; j < n; j++) {
            tmp.push_back(res[i - 1][j] + tau / 24 * (55 * f[j](time[i - 1], res[i - 1]) - 59 * f[j](time[i - 2], res[i - 2]) + 37 * f[j](time[i - 3], res[i - 3]) - 9 * f[j](time[i - 4], res[i - 4])));
        }
        res.push_back(tmp);
    }
    return res;
}



vector<vector<double>> PredictionCorection(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
    // Сетка
    int n = f.size();
    vector<double> time = {};
    double t0 = 0;
    while (t0 <= T) {
        time.push_back(t0);
        t0 += tau;
    }


    // Найдем y1, y2, y3 с помощью метода Эйлера
    vector<vector<double>> res{ y0 };

    for (int i = 1; i < 4; i++) {
        vector<double> tmp{};

        for (int j = 0; j < n; j++) {
            tmp.push_back(res[i - 1][j] + tau * f[j](time[i - 1], res[i - 1]));
        }
        res.push_back(tmp);
    }

    // Продолжим с помощью метода Адамса Бэшфорда 4-го порядка
    for (int i = 4; i < time.size(); i++) {
        // Прогноз
        vector<double> tmp{};
        for (int j = 0; j < n; j++) {
            tmp.push_back(res[i - 1][j] + tau / 24 * (55 * f[j](time[i - 1], res[i - 1]) - 59 * f[j](time[i - 2], res[i - 2]) + 37 * f[j](time[i - 3], res[i - 3]) - 9 * f[j](time[i - 4], res[i - 4])));
        }
        res.push_back(tmp);

        // Коррекция
        for (int j = 0; j < n; j++) {
            res[i][j] = res[i - 1][j] + tau / 24 * (9 * f[j](time[i], tmp) + 19 * f[j](time[i - 1], res[i - 1]) - 5 * f[j](time[i - 2], res[i - 2]) + f[j](time[i - 3], res[i - 3]));
        }
    }
    return res;
}

void ProcessAitken(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f, double q, vector<vector<double>>(&solver)(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f)) {
    vector<vector<double>> sol1 = solver(tau, T, y0, f);
    vector<vector<double>> sol2 = solver(tau * q, T, y0, f);
    vector<vector<double>> sol3 = solver(tau * q * q, T, y0, f);
    int center = sol1.size() / 10;
    cout << log(abs((sol3[center / (q * q)][0] - sol2[center / q][0]) / (sol2[center / q][0] - sol1[center][0]))) / log(q);
}

vector<vector<double>> Euler_explicit(
    vector<double> y0,
    const vector<double>& grid,
    vector<double(*)(double, vector<double>)> syst) {

    int n = grid.size(); // количество узлов
    int num_vars = syst.size(); // количество переменных
    double tau = grid[1] - grid[0];
    vector<vector<double>> res; // res хранит вектор по каждой слою времени, 
    // в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i

    res.emplace_back(y0); // добавляем начальный слой времени

    // i - номер узла, j - номер переменной
    for (int i = 1; i < n; i++) {
        vector<double> y_cur(num_vars);
        for (int j = 0; j < num_vars; j++) {
            y_cur[j] = res[i - 1][j] + tau * syst[j](grid[i-1], res[i-1]);

        }
        res.emplace_back(y_cur);
    }

    return res;
}

void WriteImplicitEuler(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
    vector<vector<double>> Matrix = ImplicitEuler(tau, T, y0, f);
    ofstream outFile("ImplicitEuler.txt");
    for (auto& row : Matrix) {
        for (auto& el : row) {
            outFile << el << " ";
        }
        outFile << endl;
    }
    outFile.close();

}


void WriteAdamsBashford(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
    vector<vector<double>> Matrix = AdamsBashford(tau, T, y0, f);
    ofstream outFile("AdamsBashford.txt");
    for (auto& row : Matrix) {
        for (auto& el : row) {
            outFile << el << " ";
        }
        outFile << endl;
    }
    outFile.close();
}


void WritePredictionCorection(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
    vector<vector<double>> Matrix = PredictionCorection(tau, T, y0, f);
    ofstream outFile("PredictionCorection.txt");
    for (auto& row : Matrix) {
        for (auto& el : row) {
            outFile << el << " ";
        }
        outFile << endl;
    }
    outFile.close();
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


vector<vector<double>> SymmetricalScheme(vector<double> y0,
    const vector<double>& grid,
    vector<double(*)(double, vector<double>)> syst) {
    int n = grid.size(); // количество узлов
    int num_vars = syst.size(); // количество переменных
    double tau = grid[1] - grid[0];
    vector<vector<double>> res; // res хранит вектор по каждой слою времени, 
    // в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i

    res.emplace_back(y0);
    for (int i = 1; i < n; i++) {
        vector<double> y_cur;
        vector<function<double(vector<double>)>> functions; // система нелинейных уравнений
        for (int j = 0; j < num_vars; j++) {
            auto func{ [tau, i, j,res,syst, grid](vector<double> coordinates)
                {
                    double t = coordinates[0]; // syst[j] принимает на вход t, y, поэтому выделяем t
                    coordinates.erase(coordinates.begin());
                    double ans = ((coordinates[j] - res[i - 1][j]) / tau);

                    ans -= 0.5 * (syst[j](grid[i - 1], res[i - 1]) + syst[j](t, coordinates));
                    return ans;
                }
            };

            functions.emplace_back(func);

        }
        y_cur = NewtonMethod(grid[i], res[i - 1], functions); // в качестве времени считаем текущее grid[i], а за начальное приближение считаем предыдущее решение res[i - 1]
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
    /// Пример пружина неявный Эйлер
    /*WriteImplicitEuler(0.1, 10, { 2,0 }, { f1_system0_test , f2_system0_test });*/

    // Пример аналитический 1
    /*WriteImplicitEuler(0.1, 1, { 2,0 }, { f1_system1_test , f2_system1_test });
    WriteAdamsBashford(0.1, 1, { 2,0 }, { f1_system1_test , f2_system1_test });
    WritePredictionCorection(0.1, 1, { 2,0 }, { f1_system1_test , f2_system1_test });*/

    // Пример аналитический 2
    /*WriteImplicitEuler(0.1, 1, { 2,0 }, { f1_system2_test , f2_system2_test });
    WriteAdamsBashford(0.1, 1, { 2,0 }, { f1_system2_test , f2_system2_test });
    WritePredictionCorection(0.1, 1, { 2,0 }, { f1_system2_test , f2_system2_test });*/

    // Порядок Эйлер меняем тау от 0.1 до 0.01 и видно стремление к 1

    /*ProcessAitken(0.1, 10, { 2,0 }, { f1_system1_test , f2_system1_test },0.5, ImplicitEuler);*/

    // Порядок AdamsBashford меняем тау от 0.1 до 0.01 и видно стремление к 1
    /*ProcessAitken(0.1, 10, { 2,0 }, { f1_system1_test , f2_system1_test }, 0.5, AdamsBashford);*/

    // Порядок PredictionCorection меняем тау от 0.1 до 0.01 и видно стремление к 1
    ProcessAitken(0.1, 10, { 2,0 }, { f1_system1_test , f2_system1_test }, 0.5, PredictionCorection);

    WriteImplicitEuler(0.001, 10, { 0,0 }, { f1_system1_book , f2_system1_book });
    vector<double> grid = GenerateUniformGrid(0, 1, 10);
    vector<vector<double>> res = Euler_explicit({ 2,0 }, grid, { f1_system1_test, f2_system1_test });
    PrintGridFunc(res);

    res = SymmetricalScheme({ 2,0 }, grid, { f1_system1_test, f2_system1_test });
    //vector<function<int(int)>> vec;
    PrintGridFunc(res);
}
