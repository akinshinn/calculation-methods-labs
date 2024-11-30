// lab5.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

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
const int BINS = 100;
const int n1 = 10;
const int n2=10;
const double sn = 1e-8;  // sn - Small Number

double test1(double x) {
    return (x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75);
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


//vector<double> NewtonMethod(double a, double b, double(&func)(double x)) {
//    pair<vector<double>, vector<double>> interval_roots = LocalizeRoots(a, b, BINS, func);
//    int n = interval_roots.first.size();
//    int iterations = 0;
//    vector<double> roots(n);
//    double x1, x2;
//    for (int i = 0; i < n; ++i) {
//        x1 = interval_roots.first[i];
//        do {
//            x2 = x1 - 
//        } while ((abs(x1 - x2) > eplsilon) && (iterations < MAX_ITER))
//    }
//    cout << "Iterations in Newton method = " << iterations << endl;
//    return roots;
//}


// varibale -- переменная по которой будет взята частная производная
double fun1(double x1, double x2) {
    return x1 * x1 - x2 * x2 - 15;
}

double fun2(double x1, double x2) {
    return x1 * x2 + 4;
}

double Derivativefun1(double x1, double x2, int var) {
    if (var == 1) {return 2 * x1; }
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

int main()
{
    cout << "test1" << endl;
    DisplayVector(BisectionMethod(0, 1, test1));
    cout << endl;
    cout << "test4" << endl;
    cout << "Derivative approximately" << endl;
    Display(SearchRoots2D(10, 10, fun1, fun2));

    cout << "Derivative accurate" << endl;
    Display(AccurateSearchRoots2D(10, 10, fun1, fun2, Derivativefun1, Derivativefun2));

    cout << endl;
    cout << "test5" << endl;
    cout << "Derivative approximately" << endl;
    Display(SearchRoots2D(10, 10, fun3, fun4));

    cout << "Derivative accurate" << endl;
    Display(AccurateSearchRoots2D(10, 10, fun3, fun4, Derivativefun3, Derivativefun4));
}
