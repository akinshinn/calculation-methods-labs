#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
using namespace std;
const int step_t = 1;
const double pi = 3.141592653589793238462643383279;

double f1(double x) {
    return sin(pi * x);
}
double f1SecondDer(double x) {
    return -pi* pi * sin(pi * x);
}

double g1(double x) {
    return 0;
}
double f2(double x) {
    return x * (1 - x);
}
double f2SecondDer(double x) {
    return -2;
}
void DisplayVector(vector<double> vec) {
    for (auto el : vec) {
        cout << el << " ";
    }
    cout << endl;
}

void WritingVectorFile(string filename, vector<double> vec) {
    ofstream file(filename, ios::app);
    for (auto& el : vec) {
        file <<setprecision(10) << el << " ";
    }
    file << endl;
}

void ExplicitMethodKolebania(string filename, double L, double T, double tau, double h, double a, double(&InitialPositionString)(double x),
    double(&InitialPositionStringSecondDerivative)(double x), double(&InitialVelocityString)(double x),
    double(&leftBoundary)(double x), double(&RightBoundary)(double x)) {
    ofstream file(filename);
    // делаем сетку
    vector<double> gridt{0}, gridx{0};
    vector<double> prev{},cur{}, next{};
    double t = 0, x=0;

    while (abs(t - T) >= 1e-10) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (abs(x-L) >= 1e-10) {
        x += h;
        gridx.emplace_back(x);
    }
    WritingVectorFile("gridx.txt", gridx);
    WritingVectorFile("gridt.txt", gridt);

    int n = gridx.size(), m = gridt.size();
    // Заполняем 0ой слой и 1ый
    for (int i = 0; i < n; i++) {
        prev.push_back(InitialPositionString(gridx[i]));
        cur.push_back(prev[i] + InitialVelocityString(gridx[i]) * tau + a * a * tau * tau * 0.5 * InitialPositionStringSecondDerivative(gridx[i]));
    }
    WritingVectorFile(filename,prev);
    if (step_t == 1) { WritingVectorFile(filename, cur); }

    DisplayVector(prev);
    for (int j = 2; j < m; j++) {
        // Левые граничные условия
        next = {};
        next.push_back(leftBoundary(gridt[j]));
        //Внутренние точки
        for (int i = 1; i < n - 1; i++) {
            next.push_back(tau * tau * a * a / (h * h) * (cur[i + 1] - 2 * cur[i] + cur[i - 1]) + 2 * cur[i] - prev[i]);
        }
        //Правые граничные 
        next.push_back(RightBoundary(gridt[j]));
        if (j % step_t==0) { WritingVectorFile(filename, next); }
        prev = cur;
        cur = next;
    }
}
int main()
{
    // Пример 1
    ExplicitMethodKolebania("Example1.txt", 5, 100, 0.01, 0.1, 1, f1, f1SecondDer, g1, g1, g1);
    /*double h = 0.1, x = 0, L = 5;
    vector<double> gridx{};
    while (abs(x - L) >= 1e-10) {
        x += h;
        gridx.emplace_back(x);
    }
    for (auto& el : gridx) {
        cout << f1(el) << endl;
    }*/

}

