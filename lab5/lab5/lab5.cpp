#include <fstream>
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
const double epsilon = 1e-6;
const int MAX_ITER = 100;
const int BINS = 73;
const int n1 = 10;
const int n2=10;
const double sn = 1e-8;  // sn - Small Number



double var1(double x) {
    double R1 = exp((x * x * x * x + x * x - x + sqrt(5)) / 5);
    double R2 = sinh((x * x * x + 21 * x + 9) / (21 * x + 6));
    return R1 + R2 - 3;
}


double der_var1(double x) {
    return 0.2 * exp(0.2* (sqrt(5) -x + x*x + x * x *x *x))* (-1 +2*x +4*x*x*x)+ 
        ((21+3*x*x)/(6+21*x) - (21*(9+21*x +x*x*x))/((6+21*x)* (6 + 21 * x)))*cosh((9+21*x+x*x*x)/(6+21*x));
}


double f1(double x, double y) {
    return 2 * x + y * y * y - 3;
}


double der_f1(double x, double y,int var) {
    if (var == 1) return 2;
    return 3 * y * y;
}


double f2(double x, double y) {
    return 2 * x * x + y * y - 3;
}


double test_newton1(double x, double y) {
    return 2 * (x + y) * (x + y) + (x - y) * (x - y) - 8;
}


double test_newton2(double x, double y) {
    return 5*x*x+(y-3)*(y-3)-9;
}

double der_f2(double x, double y, int var) {
    if (var == 2)
    return 2 * y;
    return 4 * x;
}


double x_sqr(double x) {
    return x * x;
}


double der_x_sqr(double x) {
    return 2*x;
}

double test1(double x) {
    return (x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75);
}


double exact_der_test1(double x) {
    return (-0.75 + x) * (-0.7 + x) * (-0.55 + x) * (-0.22 + x) + (-0.75 + x) * (-0.7 +
        x) * (-0.55 + x) * (-0.1 + x) + (-0.75 + x) * (-0.7 + x) * (-0.22 +
            x) * (-0.1 + x) + (-0.75 + x) * (-0.55 + x) * (-0.22 + x) * (-0.1 +
                x) + (-0.7 + x) * (-0.55 + x) * (-0.22 + x) * (-0.1 + x);
}
double test2(double x) {
    return sqrt(x+1) - 1;
}

double exact_der_test2(double x) {
    return 1 / (2 * sqrt(x + 1));
}

double test3(double x) {
    return 35 * x * x * x - 67 * x * x - 3 * x + 3;
}


double exact_der_test3(double x) {
    return -3 - 134*x + 105*x *x;
}
// varibale -- переменная по которой будет взята частная производная
double fun1(double x1, double x2) {
    return x1 * x1 - x2 * x2 - 15;
}

double fun2(double x1, double x2) {
    return x1 * x2 + 4;
}

double Derivativefun1(double x1, double x2, int var) {
    if (var == 1) { return 2 * x1; }
    return -2 * x2;
}

double Derivativefun2(double x1, double x2, int var) {
    if (var == 1) { return x1; }
    return x2;
}

double fun3(double x1, double x2) {
    return x1 * x1 + x2 * x2 + x1 + x2 - 8;
}

double fun4(double x1, double x2) {
    return x1 * x1 + x2 * x2 + x1 * x2 - 7;
}

double Derivativefun3(double x1, double x2, int var) {
    if (var == 1) { return 2 * x1 + 1; }
    return 2 * x2 + 1;
}

double Derivativefun4(double x1, double x2, int var) {
    if (var == 1) { return 2 * x1 + x2; }
    return 2 * x2 + x1;
}

double Derivative(pair<double, double> Point, double(&func)(double x1, double x2), int variable) {
    if (variable == 1) {
        return (func(Point.first + sn, Point.second) - func(Point.first, Point.second)) / sn;
    }
    return (func(Point.first, Point.second + sn) - func(Point.first, Point.second)) / sn;

}

double Derivative1D(double x, double(&func)(double x)) {
    //return (func(Point.first, Point.second + sn) - func(Point.first, Point.second)) / sn;
    return (func(x + sn) - func(x - sn)) / (2 * sn);

}

double Derivative1D(double x, double FixedArg, int num, double(&func)(double x, double y)) {
    //return (func(Point.first, Point.second + sn) - func(Point.first, Point.second)) / sn;
    if (num == 1) {
        return (func(x + sn, FixedArg) - func(x - sn, FixedArg)) / (2 * sn);
    }
    return (func(FixedArg, x + sn) - func(FixedArg, x - sn)) / (2 * sn);

}
void DisplayMatrix(vector<vector<double>> Matrix) {
    for (auto& row : Matrix) {
        for (auto el : row) {
            cout << el << " ";
        }
        cout << endl;
    }
}

void Display(vector<pair<double, double>> vec) {
    for (auto& el : vec) {
        cout << el.first << " " << el.second << endl;
    }
}


void AddElement(vector<pair<double, double>>& vec, pair<double, double> adpar) {
    int flag = 0;
    
    for (auto& el : vec) {
        if ((el.first - adpar.first) < 1e-4 && (el.second-adpar.second) < 1e-4) { flag++; break; }
    }
    if (flag == 0) {
        vec.push_back(adpar);
    }
}


void DisplayVector(vector<double> vec) {
    for (auto el : vec) {
        cout << el << " ";
    }
    cout << endl;
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


pair<vector<double>, vector<double>> GenerateUniformGrid(double a, double b, int n, double(&func)(double x)) {
    // n - количество узлов, [a, b] - отрезок интерполирования
    vector<double> xValues{}, yValues{};
    double h = (b - a) / n;
    for (int i = 0; i < n + 1; i++) {
        xValues.push_back(a + i * h);
        yValues.push_back(func(a + i * h));
    }
    return { xValues, yValues };
}

pair<vector<pair<double, double>>, vector<pair<double, double>>> LocalizeRoots(double a, double b, int size, double(&func)(double x)) {
    pair<vector<double>, vector<double>> grid = GenerateUniformGrid(a, b, size, func);
    vector<double> x = grid.first;
    vector<double> y = grid.second;
    int n = grid.first.size();
    vector<pair<double, double>> left_boundry, right_boundry;
    int i = 0;
    while (i < n-1) {
        if (y[i] * y[i + 1] <= 0) {
            left_boundry.push_back({ x[i],y[i]});
            right_boundry.push_back({x[i+1],y[i+1]});
        }
        i++;
    }
    return { left_boundry, right_boundry };
}


vector<double> BisectionMethod(double a, double b, double(&func)(double x)) {
    pair<vector<pair<double, double>>, vector<pair<double, double>>> interval_roots = LocalizeRoots(a, b, BINS, func);
    int n = interval_roots.first.size();
    int iterations = 0;
    vector<double> roots(n);
    for (int i = 0; i < n; ++i) {
        double x1, x2, y1, y2, mid, ymid;
        x1 = interval_roots.first[i].first;
        y1 = interval_roots.first[i].second;
        x2 = interval_roots.second[i].first;
        y2 = interval_roots.second[i].second;

        if (y1 == 0) {
            roots[i] = x1;
        }
        else if (y2 == 0) {
            roots[i] = x2;
        }
        else {
            while (((x2 - x1) > 2*epsilon) && (iterations < MAX_ITER)) {
                mid = (x1 + x2) / 2;
                ymid = func(mid);
                if ((y1 < 0) && (y2 > 0)) {
                    if (ymid == 0) { roots[i] = mid; break; }
                    else if (ymid > 0) { x2 = mid; y2 = ymid; }
                    else { x1 = mid; y1 = ymid; }
                }
                else { 
                    if (ymid == 0) { roots[i] = mid; break; }
                    else if (ymid > 0) { x1 = mid; y1 = ymid; }
                    else { x2 = mid; y2 = ymid; }
                }
                iterations++;
            }
            roots[i] = (x1 + x2) / 2;
        }

    }
    cout << "Iterations in bisection method = " << iterations << endl;
    //фильтрация корней тк если корни попали на границу сетки то один корень учитавается 2 раза.
    int i=0, cntRoot = 0;
    vector<double> Roots = {};
    while (i < n - 1) {
        if (roots[i] == roots[i + 1]) {
            cntRoot++;
            Roots.push_back(roots[i]);
            i += 2;
        }
        else {
            cntRoot++;
            Roots.push_back(roots[i]);
            i++;
        }
    }
    if (i <= n - 1) {
        Roots.push_back(roots[i]);
    }
    return Roots;
}


// Общий случай метода Ньютона с n корнями и численной производной
vector<double> NewtonMethod(double a, double b, double(&func)(double x)) {
    pair<vector<pair<double, double>>, vector<pair<double, double>>> interval_roots = LocalizeRoots(a, b, BINS, func);
    int n = interval_roots.first.size();
    int iterations = 0;
    vector<double> roots(n);
    double x1, x2, y1, y2;
    for (int i = 0; i < n; ++i) {
        y1 = interval_roots.first[i].second;
        x2 = interval_roots.first[i].first; // Чтобы в начале цикла x1 = x2, все равно потом x2 поменяется
        y2 = interval_roots.second[i].second;
        if (y1 == 0) { roots[i] = interval_roots.first[i].first; continue; }
        if (y2 == 0) { roots[i] = interval_roots.second[i].first; continue; }
        do {
            iterations++;
            x1 = x2;
            double der_value = Derivative1D(x1, func);
            x2 = x1 - (func(x1)) / (der_value);
        } while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));
        roots[i] = x2;
    }
    int i = 0, cntRoot = 0;
    vector<double> Roots = {};
    while (i < n - 1) {
        if (roots[i] == roots[i + 1]) {
            cntRoot++;
            Roots.push_back(roots[i]);
            i += 2;
        }
        else {
            cntRoot++;
            Roots.push_back(roots[i]);
            i++;
        }
    }
    if (i <= n - 1) {
        Roots.push_back(roots[i]);
    }

    cout << "Iterations in Newton method = " << iterations << endl;
    return Roots;
}


// Общий случай метода Ньютона с n корнями и точной производной
vector<double> NewtonMethod(double a, double b, double(&func)(double x), double(&der)(double x)) {
    pair<vector<pair<double, double>>, vector<pair<double, double>>> interval_roots = LocalizeRoots(a, b, BINS, func);
    int n = interval_roots.first.size();
    int iterations = 0;
    vector<double> roots(n);
    double x1, x2, y1, y2;
    for (int i = 0; i < n; ++i) {
        y1 = interval_roots.first[i].second;
        x2 = interval_roots.first[i].first; // Чтобы в начале цикла x1 = x2, все равно потом x2 поменяется
        y2 = interval_roots.second[i].second;
        if (y1 == 0) { roots[i] = interval_roots.first[i].first; continue; }
        if (y2 == 0) { roots[i] = interval_roots.second[i].first; continue; }
        do {
            iterations++;
            x1 = x2;
            double der_value = der(x1);
            x2 = x1 - (func(x1)) / (der_value);
        } while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));
        roots[i] = x2;
    }
    int i = 0, cntRoot = 0;
    vector<double> Roots = {};
    while (i < n - 1) {
        if (roots[i] == roots[i + 1]) {
            cntRoot++;
            Roots.push_back(roots[i]);
            i += 2;
        }
        else {
            cntRoot++;
            Roots.push_back(roots[i]);
            i++;
        }
    }
    if (i <= n - 1) {
        Roots.push_back(roots[i]);
    }

    cout << "Iterations in Newton method = " << iterations << endl;
    return Roots;
}


// Случай с 1 корнем и численной производной
vector<double> NewtonMethod(double a, double b, double x0, double(&func)(double x)) {
    int n = 1;
    int iterations = 0;
    vector<double> roots(n);
    double x1 = x0, x2 = x0;
    double lastx;
    for (int i = 0; i < n; ++i) {
        do {
            iterations++;
            lastx = x1;
            x1 = x2;
            double der_value = Derivative1D(x1, func);
            x2 = x1 - (func(x1)) / (der_value);
            cout << iterations << " & " << abs(x2 - x1) << " & " << abs(x2 - 0.485204) << " & " <<
                log(abs((x2 - 0.485204) / (x1 - 0.485204))) / log(abs((x1 - 0.485204) / (lastx - 0.485204))) << " \\\\ \\hline" << endl;
            /*} while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));*/ // Оригинальное условие
        } while (iterations < 20); // условие для оценки скорости
        roots[i] = x2;
    }

    cout << "Iterations in Newton method = " << iterations << endl;
    return roots;
}


// Случай с 1 корнем и точной производной
vector<double> NewtonMethod(double a, double b, double x0, double(&func)(double x), double(&der)(double x)) {
    int n = 1;
    int iterations = 0;
    vector<double> roots(n);
    double x1 = x0 , x2 = x0;
    double lastx;
    for (int i = 0; i < n; ++i) {
        do {
            iterations++;
            lastx = x1;
            x1 = x2;
            double der_value = der(x1);
            x2 = x1 - (func(x1)) / (der_value);
            //cout << iterations << " & " << abs(x2 - x1) << " & " << abs(x2 - 0.485204) << " & " <<
            //    log(abs((x2 - 0.485204) / (x1 - 0.485204))) / log(abs((x1 - 0.485204) / (lastx - 0.485204))) << " \\\\ \\hline" << endl;
            /*} while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));*/ // Оригинальное условие
        } while (iterations < 20); // условие для оценки скорости
        roots[i] = x2;
    }

    cout << "Iterations in Newton method = " << iterations << endl;
    return roots;
}


// Модификация (исключение выхода за границы отрезка) метода Ньютона с 1 корнем и численной производной
vector<double> NewtonMethodModification(double a, double b, double x0, double(&func)(double x)) {
    int n = 1;
    int iterations = 0;
    vector<double> roots(n);
    double x1 = x0, x2 = x0;
    double lastx;
    for (int i = 0; i < n; ++i) {
        do {
            iterations++;
            lastx = x1;
            x1 = x2;
            double der_value = Derivative1D(x1, func);
            x2 = x1 - (func(x1)) / (der_value);
            if (x2 < a) x2 = (a * func(b) - b*func(a))/(func(b) - func(a)) ;
            else if (x2 > b) x2 = (a * func(b) - b * func(a)) / (func(b) - func(a));
            cout  << iterations << " & " << abs(x2 - x1) << " & " << abs(x2 - 0.485204) << " & " << log(abs((x2 - x1) / (x1 - 0.485204))) / log(abs((x1 - 0.485204) / (x1 - lastx))) << "\\ \hline" << endl;
            /*} while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));*/ // Оригинальное условие
        } while (iterations < 20); // условие для оценки скорости
        roots[i] = x2;
    }

    cout << "Iterations in Newton method = " << iterations << endl;
    return roots;
}


// Модификация (исключение выхода за границы отрезка и зацикливание) метода Ньютона с 1 корнем и точной производной
vector<double> NewtonMethodModification(double a, double b, double x0, double(&func)(double x), double(&der)(double x)) {
    int n = 1;
    int iterations = 0;
    vector<double> roots(n);
    double x1, x2 = x0;
    for (int i = 0; i < n; ++i) {
        do {
            iterations++;
            x1 = x2;
            double der_value = der(x1);
            //cout << "der = " << der_value << endl;
            x2 = x1 - (func(x1)) / (der_value)+0.5 * epsilon;
            if (x2 < a) x2 = (a * func(b) - b * func(a)) / (func(b) - func(a));
            else if (x2 > b) x2 = (a * func(b) - b * func(a)) / (func(b) - func(a));
            //cout << "x2 = " << x2 << endl;
        } while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));
        roots[i] = x2;
    }

    cout << "Iterations in Newton method = " << iterations << endl;
    return roots;
}

vector<double> NewtonMethodModificationForSystem(double a, double b, double x0,double FixedVar, int num,  double(&func)(double x, double y)) {
    int n = 1;
    int iterations = 0;
    vector<double> roots(n);
    double x1 = x0, x2 = x0;
    double lastx;
    for (int i = 0; i < n; ++i) {
        do {
            iterations++;
            lastx = x1;
            x1 = x2;
            double der_value = Derivative1D(x1, FixedVar, num, func);
            if (num == 1) {
                x2 = x1 - (func(x1, FixedVar)) / (der_value);
                if (x2 < a) x2 = (a * func(b, FixedVar) - b * func(a, FixedVar)) / (func(b, FixedVar) - func(a, FixedVar));
                else if (x2 > b) x2 = (a * func(b, FixedVar) - b * func(a, FixedVar)) / (func(b, FixedVar) - func(a, FixedVar));
                //cout << iterations << " & " << abs(x2 - x1) << " & " << abs(x2 - 0.485204) << " & " << log(abs((x2 - x1) / (x1 - 0.485204))) / log(abs((x1 - 0.485204) / (x1 - lastx))) << "\\ \hline" << endl;
                /*} while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));*/ // Оригинальное условие
            }
            else {
                x2 = x1 - (func(FixedVar,x1)) / (der_value);
                if (x2 < a) x2 = (a * func(FixedVar, b) - b * func(FixedVar, a)) / (func(FixedVar, b) - func(FixedVar, a));
                else if (x2 > b) x2 = (a * func(FixedVar, b) - b * func(FixedVar, a)) / (func(FixedVar, b) - func(FixedVar, a));
                //cout << iterations << " & " << abs(x2 - x1) << " & " << abs(x2 - 0.485204) << " & " << log(abs((x2 - x1) / (x1 - 0.485204))) / log(abs((x1 - 0.485204) / (x1 - lastx))) << "\\ \hline" << endl;
            }
        } while (iterations < 20); // условие для оценки скорости
        roots[i] = x2;
    }

    //cout << "Iterations in Newton method = " << iterations << endl;
    return roots;
}

vector<double> NewtonMethodModificationForSystem(double a, double b, double x0,double FixedVar,int num, double(&func)(double x, double y), double(&der)(double x, double y, int var)) {
    int n = 1;
    int iterations = 0;
    vector<double> roots(n);
    double x1, x2 = x0;
    for (int i = 0; i < n; ++i) {
        do {
            iterations++;
            x1 = x2;
            if (num==1){
                double der_value = der(x1, FixedVar, 1);
                x2 = x1 - (func(x1, FixedVar)) / (der_value)+0.5 * epsilon;
                if (x2 < a) x2 = (a * func(b, FixedVar) - b * func(a, FixedVar)) / (func(b, FixedVar) - func(a, FixedVar));
                else if (x2 > b) x2 = (a * func(b, FixedVar) - b * func(a, FixedVar)) / (func(b, FixedVar) - func(a, FixedVar));
            }
            else {
                double der_value = der(FixedVar, x1, 2);
                x2 = x1 - (func(FixedVar, x1)) / (der_value)+0.5 * epsilon;
                if (x2 < a) x2 = (a * func(FixedVar, b) - b * func(FixedVar, a)) / (func(FixedVar, b) - func(FixedVar, a));
                else if (x2 > b) x2 = (a * func(FixedVar, b) - b * func(FixedVar, a)) / (func(FixedVar, b) - func(FixedVar, a));
            }
            //cout << "der = " << der_value << endl;

            //cout << "x2 = " << x2 << endl;
        } while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));
        roots[i] = x2;
    }

    //cout << "Iterations in Newton method = " << iterations << endl;
    return roots;
}

vector<double> SecantMethod(double a, double b, double(&func)(double x)) {
    pair<vector<pair<double, double>>, vector<pair<double, double>>> interval_roots = LocalizeRoots(a, b, BINS, func);
    int n = interval_roots.first.size();
    int iterations = 0;
    vector<double> roots(n);
    double x1, x2, y1, y2;
    for (int i = 0; i < n; ++i) {
        y1 = interval_roots.first[i].second;
        x2 = interval_roots.first[i].first; // Чтобы в начале цикла x1 = x2, все равно потом x2 поменяется
        y2 = interval_roots.second[i].second;
        if (y1 == 0) { roots[i] = interval_roots.first[i].first; continue; }
        if (y2 == 0) { roots[i] = interval_roots.second[i].first; continue; }
        b = interval_roots.second[i].first;
        a = interval_roots.first[i].first;
        do {
            iterations++;
            x1 = x2;
            
            x2 = x1 - func(x1)* (b-x1) / (func(b) - func(x1));
        } while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));
        roots[i] = x2;
    }
    int i = 0, cntRoot = 0;
    vector<double> Roots = {};
    while (i < n - 1) {
        if (roots[i] == roots[i + 1]) {
            cntRoot++;
            Roots.push_back(roots[i]);
            i += 2;
        }
        else {
            cntRoot++;
            Roots.push_back(roots[i]);
            i++;
        }
    }
    if (i <= n - 1) {
        Roots.push_back(roots[i]);
    }

    cout << "Iterations in Secant method = " << iterations << endl;
    return Roots;
}


double SecantMethodOneRoot(double a, double b, double x0, double(&func)(double x)) {
    int iterations = 0;
    double x1 = (b-a)/2 , x2 = (b - a) / 2;
    double lastx=a;
    double res;
    do {
        iterations++;


        x2 = lastx - func(lastx) * (x1 - lastx) / (func(x1) - func(lastx));
        cout << iterations << " & " << abs(x2 - x1) << " & " << abs(x2 - 0.485204) << " & " << 
            log(abs((x2 - 0.485204) / (x1 - 0.485204))) / log(abs((x1 - 0.485204) / (lastx - 0.485204))) << " \\\\ \\hline" << endl;
        //cout << "x2 = " << x2 << endl;
        lastx = x1;
        x1 = x2;
        /*} while ((abs(x1 - x2) > epsilon) && (iterations < MAX_ITER));*/ // Оригинальное условие
    } while (iterations < 20); // условие для оценки скорости

    res = x2;


    cout << "Iterations in Secant method = " << iterations << endl;
    return res;
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

// возвращает два приближения последних для того тчобы понять алгоритм сошелся или по количество закончили
pair<pair<double, double>, pair<double, double>>NewtonMethodSystem(pair<double, double> StartPoint, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2), int& iteration) {
    int iterations = 0;
    pair<double, double> nextPoint= StartPoint, curPoint = StartPoint;
    vector<double> sol{};
    while (iterations < MAX_ITER)
    {
        curPoint = nextPoint;
        vector<vector<double>> Jacoby{ {Derivative(curPoint, func1, 1), Derivative(curPoint, func1, 2), -func1(curPoint.first, curPoint.second)},{Derivative(curPoint, func2, 1), Derivative(curPoint, func2, 2), -func2(curPoint.first, curPoint.second)} };
        sol = Gauss_method(Jacoby);
        nextPoint.first = curPoint.first + sol[0];
        nextPoint.second = curPoint.second + sol[1];
        if (sqrt(sol[0] * sol[0] + sol[1] * sol[1]) < epsilon) { break; }
        iteration++;
        iterations++;
    }
    return { nextPoint, curPoint };
}

pair<double, double> NewtonMethodSystemModification(pair<double, double> boundar,pair<double, double> StartPoint, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2), int& iteration) {
    int iterations = 0;
    pair<double, double> nextPoint = StartPoint, curPoint = StartPoint;
    vector<double> sol{};
    while (iterations < MAX_ITER)
    {
        curPoint = nextPoint;
        vector<double> res = NewtonMethodModificationForSystem(boundar.first, boundar.second, curPoint.first, curPoint.second, 1, func1);
        nextPoint.first = res[0];
        res = NewtonMethodModificationForSystem(boundar.first, boundar.second, nextPoint.first, nextPoint.second, 2, func2);
        nextPoint.second = res[0];
        if (sqrt((nextPoint.first - curPoint.first) * (nextPoint.first - curPoint.first) + (nextPoint.second - curPoint.second) * (nextPoint.second - curPoint.second)) < epsilon) { break; }
        iteration++;
        iterations++;
    }
    return { nextPoint };
}

pair<double, double> NewtonMethodSystemModification(pair<double, double> boundar, pair<double, double> StartPoint, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2), double(&der1)(double x1, double x2, int var), double(&der2)(double x1, double x2, int var), int& iteration) {
    int iterations = 0;
    pair<double, double> nextPoint = StartPoint, curPoint = StartPoint;
    vector<double> sol{};
    while (iterations < MAX_ITER)
    {
        curPoint = nextPoint;
        vector<double> res = NewtonMethodModificationForSystem(boundar.first, boundar.second, curPoint.first, curPoint.second, 1, func1 , der1);
        nextPoint.first = res[0];
        res = NewtonMethodModificationForSystem(boundar.first, boundar.second, nextPoint.first, nextPoint.second, 2, func2, der2);
        nextPoint.second = res[0];
        if (sqrt((nextPoint.first - curPoint.first) * (nextPoint.first - curPoint.first) + (nextPoint.second - curPoint.second) * (nextPoint.second - curPoint.second)) < epsilon) { break; }
        iteration++;
        iterations++;
    }
    return { nextPoint };
}
///////////////////////////////////////


pair<pair<double, double>, pair<double, double>>NewtonMethodSystemAccurate(pair<double, double> StartPoint, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2), double(&derfunc1)(double x1, double x2,int var), double(&derfunc2)(double x1, double x2,int var), int& iteration) {
    int iterations = 0;
    pair<double, double> nextPoint = StartPoint, curPoint = StartPoint;
    vector<double> sol{};
    while (iterations < MAX_ITER)
    {
        curPoint = nextPoint;
        vector<vector<double>> Jacoby{ {derfunc1(curPoint.first,curPoint.second,1),derfunc1(curPoint.first,curPoint.second,2), -func1(curPoint.first, curPoint.second)},{derfunc2(curPoint.first,curPoint.second,1), derfunc2(curPoint.first,curPoint.second,2), -func2(curPoint.first, curPoint.second)}};
        sol = Gauss_method(Jacoby);
        nextPoint.first = curPoint.first + sol[0];
        nextPoint.second = curPoint.second + sol[1];
        if (sqrt(sol[0] * sol[0] + sol[1] * sol[1]) < epsilon) { break; }
        iteration++;
        iterations++;
    }
    return { nextPoint, curPoint };
}
// n1 - колво узлов по 1ой переменной, n2 - колво узлов по 2ой переменной


vector<pair<double, double>> GenerateUniformGrid2D(double L1, double L2, int n1, int n2) {
    vector<pair<double, double>> grid;
    vector<double> tmp{};
    double h1 = (2 * L1) / n1, h2 = (2 * L2) / n2;
    for (int i = 0; i < n1 + 1; i++) {
        for (int j = 0; j < n2 + 1; j++) {
            grid.push_back({ -L1 + i * h1, -L2 * j * h2 });
        }
    }
    return grid;
}


vector<pair<double, double>> SearchRoots2D(double L1, double L2, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2)) {
    vector<pair<double, double>> grid = GenerateUniformGrid2D(L1, L2, n1, n2);
    vector<pair<double, double>> Roots{};
    int iteration = 0;
    pair<pair<double, double>, pair<double, double>> root = {};
    for (int i = 0; i < grid.size(); i++) {
        root = NewtonMethodSystem(grid[i], func1, func2, iteration);
        if (sqrt((root.first.first - root.second.first) * (root.first.first - root.second.first) + (root.first.second - root.second.second) * (root.first.second - root.second.second)) < epsilon) { AddElement(Roots,root.first); }
    }
    cout << "NewtonMethodSystem need " << iteration << " iterations" << endl;
    return Roots;
}

vector<pair<double, double>> AccurateSearchRoots2D(double L1, double L2, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2), double(&derfunc1)(double x1, double x2, int var), double(&derfunc2)(double x1, double x2, int var)) {
    vector<pair<double, double>> grid = GenerateUniformGrid2D(L1, L2, n1, n2);
    vector<pair<double, double>> Roots{};
    int iteration = 0;
    pair<pair<double, double>, pair<double, double>> root = {};
    for (int i = 0; i < grid.size(); i++) {
        root = NewtonMethodSystemAccurate(grid[i], func1, func2,derfunc1,derfunc2, iteration);
        if (sqrt((root.first.first - root.second.first) * (root.first.first - root.second.first) + (root.first.second - root.second.second) * (root.first.second - root.second.second)) < epsilon) { AddElement(Roots, root.first); }
    }
    cout << "NewtonMethodSystem need " << iteration << " iterations" << endl;
    return Roots;
}


vector<pair<double, double>> generate_2d_grid(double x0, double x1, double y0, double y1, int nx, int ny) {
    double hx = abs(x1 - x0) / nx, hy = abs(y1 - y0) / ny;
    double cx , cy;
    vector<pair<double, double>> res;
    for (int i = 0; i <= ny; i++) {
        cx = x0;
        for (int j = 0; j <= nx; j++) {
            res.push_back({ cx,y0 });
            cout << cx << " " << y0 << endl;

            cx += hx;

        }
        y0 += hy;

    }
    return res;
}


pair<vector<pair<double, double>>, vector<int>> get_matrix_convergence(double x0, double x1, double y0, double y1, int nx, int ny, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2)) {
    vector<pair<double, double>> grid = generate_2d_grid(x0, x1, y0, y1, nx, ny);
    vector<int> all_iterations;
    for (int i = 0; i < grid.size(); ++i) {
        pair<double, double> point = grid[i];
        int cur_iter = 0;
        NewtonMethodSystem(point, func1, func2, cur_iter);
        all_iterations.push_back(cur_iter);
        cout << grid[i].first << " " << grid[i].second << " " << cur_iter << endl;
        //cout << cur_iter << endl;

    }
    return { grid, all_iterations };
}


void print_matix_convergence(pair<vector<pair<double, double>>, vector<int>> grid_iter, string file,int nx, int ny) {
    ofstream f;
    f.open(file);
    cout << grid_iter.first.size();
    if (f.is_open()) {
        for (int i = 0; i < ny; ++i) {
            for (int j = 0; j < nx; j++) {
                f << grid_iter.second[j + i] << " ";
            }
            f << endl;

        }
        f.close();
    }
}

void SearchVelocityOfConverge(pair<double, double> StartPoint, pair<double, double> answer, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2)) {
    int iterations = 0;
    pair<double, double> nextPoint = StartPoint, curPoint = StartPoint;
    vector<double> sol{};
    vector<pair<double, double>> ListSol{ curPoint,{0,0},{0,0} };
    for (int i = 1; i < 3; i++) {
        curPoint = nextPoint;
        vector<vector<double>> Jacoby{ {Derivative(curPoint, func1, 1), Derivative(curPoint, func1, 2), -func1(curPoint.first, curPoint.second)},{Derivative(curPoint, func2, 1), Derivative(curPoint, func2, 2), -func2(curPoint.first, curPoint.second)} };
        sol = Gauss_method(Jacoby);
        nextPoint.first = curPoint.first + sol[0];
        nextPoint.second = curPoint.second + sol[1];
        ListSol[i] = nextPoint;
    }
    double  ch1, ch2, ch3;
    ch1 = sqrt((ListSol[2].first - answer.first) * (ListSol[2].first - answer.first) + (ListSol[2].second - answer.second) * (ListSol[2].second - answer.second));
    ch2 = sqrt((ListSol[1].first - answer.first) * (ListSol[1].first - answer.first) + (ListSol[1].second - answer.second) * (ListSol[1].second - answer.second));
    ch3 = sqrt((ListSol[0].first - answer.first) * (ListSol[0].first - answer.first) + (ListSol[0].second - answer.second) * (ListSol[0].second - answer.second));
    cout << "k=3| " << "p = " << log(abs(ch1 / ch2)) / log(abs(ch2 / ch3)) << endl;
    while (iterations < MAX_ITER)
    {
        curPoint = nextPoint;
        vector<vector<double>> Jacoby{ {Derivative(curPoint, func1, 1), Derivative(curPoint, func1, 2), -func1(curPoint.first, curPoint.second)},{Derivative(curPoint, func2, 1), Derivative(curPoint, func2, 2), -func2(curPoint.first, curPoint.second)} };
        sol = Gauss_method(Jacoby);
        nextPoint.first = curPoint.first + sol[0];
        nextPoint.second = curPoint.second + sol[1];

        ListSol[0] = ListSol[1];
        ListSol[1] = ListSol[2];
        ListSol[2] = nextPoint;
        ch1 = sqrt((ListSol[2].first - answer.first) * (ListSol[2].first - answer.first) + (ListSol[2].second - answer.second) * (ListSol[2].second - answer.second));
        ch2 = sqrt((ListSol[1].first - answer.first) * (ListSol[1].first - answer.first) + (ListSol[1].second - answer.second) * (ListSol[1].second - answer.second));
        ch3 = sqrt((ListSol[0].first - answer.first) * (ListSol[0].first - answer.first) + (ListSol[0].second - answer.second) * (ListSol[0].second - answer.second));
        cout << "k= " << 4 + iterations << "| " << "p = " << log(abs(ch1 / ch2)) / log(abs(ch2 / ch3)) << endl;
        if (sqrt(sol[0] * sol[0] + sol[1] * sol[1]) < epsilon) { break; }
        iterations++;
    }


}


void SearchVelocityOfConvergeAccurate(pair<double, double> StartPoint, pair<double, double> answer, double(&func1)(double x1, double x2), double(&func2)(double x1, double x2), double(&derfunc1)(double x1, double x2, int var), double(&derfunc2)(double x1, double x2, int var)) {
    int iterations = 0;
    pair<double, double> nextPoint = StartPoint, curPoint = StartPoint;
    vector<double> sol{};
    vector<pair<double, double>> ListSol{ curPoint,{0,0},{0,0} };
    for (int i = 1; i < 3; i++) {
        curPoint = nextPoint;
        vector<vector<double>> Jacoby{ {derfunc1(curPoint.first,curPoint.second,1),derfunc1(curPoint.first,curPoint.second,2), -func1(curPoint.first, curPoint.second)},{derfunc2(curPoint.first,curPoint.second,1), derfunc2(curPoint.first,curPoint.second,2), -func2(curPoint.first, curPoint.second)} };
        sol = Gauss_method(Jacoby);
        nextPoint.first = curPoint.first + sol[0];
        nextPoint.second = curPoint.second + sol[1];
        ListSol[i] = nextPoint;
    }
    double  ch1, ch2, ch3;
    ch1 = sqrt((ListSol[2].first - answer.first) * (ListSol[2].first - answer.first) + (ListSol[2].second - answer.second) * (ListSol[2].second - answer.second));
    ch2 = sqrt((ListSol[1].first - answer.first) * (ListSol[1].first - answer.first) + (ListSol[1].second - answer.second) * (ListSol[1].second - answer.second));
    ch3 = sqrt((ListSol[0].first - answer.first) * (ListSol[0].first - answer.first) + (ListSol[0].second - answer.second) * (ListSol[0].second - answer.second));
    cout << "k=3| " << "p = " << log(abs(ch1 / ch2)) / log(abs(ch2 / ch3)) << endl;
    while (iterations < MAX_ITER)
    {
        curPoint = nextPoint;
        vector<vector<double>> Jacoby{ {derfunc1(curPoint.first,curPoint.second,1),derfunc1(curPoint.first,curPoint.second,2), -func1(curPoint.first, curPoint.second)},{derfunc2(curPoint.first,curPoint.second,1), derfunc2(curPoint.first,curPoint.second,2), -func2(curPoint.first, curPoint.second)} };
        sol = Gauss_method(Jacoby);
        nextPoint.first = curPoint.first + sol[0];
        nextPoint.second = curPoint.second + sol[1];

        ListSol[0] = ListSol[1];
        ListSol[1] = ListSol[2];
        ListSol[2] = nextPoint;
        ch1 = sqrt((ListSol[2].first - answer.first) * (ListSol[2].first - answer.first) + (ListSol[2].second - answer.second) * (ListSol[2].second - answer.second));
        ch2 = sqrt((ListSol[1].first - answer.first) * (ListSol[1].first - answer.first) + (ListSol[1].second - answer.second) * (ListSol[1].second - answer.second));
        ch3 = sqrt((ListSol[0].first - answer.first) * (ListSol[0].first - answer.first) + (ListSol[0].second - answer.second) * (ListSol[0].second - answer.second));
        cout << "k= " << 4 + iterations << "| " << "p = " << log(abs(ch1 / ch2)) / log(abs(ch2 / ch3)) << endl;
        if (sqrt(sol[0] * sol[0] + sol[1] * sol[1]) < epsilon) { break; }
        iterations++;
    }


}




int main()
{

    //print_matix_convergence(get_matrix_convergence(-10, 10, -10, 10, 10, 10, fun1, fun2), "matrix_conv1test.txt",10,10);
    //print_matix_convergence(get_matrix_convergence(-5, 5, -5, 5, 300, 300, test_newton1, test_newton2), "matrix_conv3test.txt", 300, 300);
    int iter = 0;
    NewtonMethodSystem({ -4,4 }, fun1, fun2,iter);
    cout << iter << endl;
    //iter = 0;
    //NewtonMethodSystem({ -10 + 2/30 * 0.5,-10 + 2/30 * 0.5 }, fun1, fun2, iter);
    //cout << iter << endl;
    //Display(SearchRoots2D(100, 100, test_newton1, test_newton2));
    //DisplayVector(NewtonMethod(0, 1, 0, var1));
    
    //cout << "test1" << endl;
    //DisplayVector(BisectionMethod(0, 1, test3));
    //cout << endl;
    //cout << "test4" << endl;
    //cout << "Derivative approximately" << endl;
    //Display(SearchRoots2D(10, 10, fun1, fun2));

    //cout << "Derivative accurate" << endl;
    //Display(AccurateSearchRoots2D(10, 10, fun1, fun2, Derivativefun1, Derivativefun2));
    //Display(AccurateSearchRoots2D(2, 2, f1, f2, der_f1, der_f2));
    //cout << endl;
    //cout << "test5" << endl;
    //cout << "Derivative approximately" << endl;
    //Display(SearchRoots2D(10, 10, fun3, fun4));

    //cout << "Derivative accurate" << endl;
    //Display(AccurateSearchRoots2D(10, 10, fun3, fun4, Derivativefun3, Derivativefun4));
    //int iter = 0;
    //pair<double, double> res = NewtonMethodSystemModification({ -10,10 }, { 3,2 }, fun1, fun2,Derivativefun1,Derivativefun2, iter);
    //cout << "sol: " << res.first << " " << res.second << ", iter =  " << iter << endl;
    //cout << "System 1" << endl;
    //SearchVelocityOfConverge({ 1,0 }, { 4,-1 }, fun1, fun2);
    //cout << endl;
    //cout << "System 2" << endl;
    //SearchVelocityOfConverge({ -2.5,0 }, { -3,1 }, fun3, fun4);
}
