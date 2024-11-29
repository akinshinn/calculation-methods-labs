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
using namespace std; 
const double eplsilon = 1e-6;
const int MAX_ITER = 100;
const int BINS = 100;


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

pair<vector<double>, vector<double>> LocalizeRoots(double a, double b, int size, double(&func)(double x)) {
    pair<vector<double>, vector<double>> grid = GenerateUniformGrid(a, b, size, func);
    vector<double> x = grid.first;
    vector<double> y = grid.second;
    int n = grid.first.size()-1;
    vector<double> left_boundry;
    vector<double> right_boundry;
    int i = 0;
    while (i < n-1) {
        if (y[i] * y[i + 1] <= 0) {
            left_boundry.emplace_back(x[i]);
            right_boundry.emplace_back(x[i+1]);
            i += 2;
            continue;
        }
        i++;
    }
    return { left_boundry, right_boundry };
}


vector<double> BisectionMethod(double a, double b, double(&func)(double x)) {
    pair<vector<double>, vector<double>> interval_roots = LocalizeRoots(a, b, BINS, func);
    int n = interval_roots.first.size();
    int iterations = 0;
    vector<double> roots(n);
    for (int i = 0; i < n; ++i) {
        double x1, x2;
        x1 = interval_roots.first[i];
        x2 = interval_roots.second[i];
        while (((x2 - x1) > eplsilon) && (iterations < MAX_ITER)) {

            x2 = (x2  + x1) / 2;
            iterations++;
        }
        roots[i] = x2;
    }
    cout << "Iterations in bisection method = " << iterations << endl;
    return roots;
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



int main()
{
    DisplayVector(BisectionMethod(0, 1, test1));
}
