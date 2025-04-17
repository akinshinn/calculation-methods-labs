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
const int tz = 1000;
double f1(double x) {
    return sin(pi * x);
}
double f1SecondDer(double x) {
    return -pi* pi * sin(pi * x);
}
double analsolEx1(double x, double t) {
    return sin(pi * x) * cos(pi * t);
}


double SecondDerNum(double(&f)(double x), double point) {
    double step = 1e-1;
    return (f(point + step) + f(point - step)-2*f(point)) / (step * step);
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
double f3(double x) {
    return x* (x + 1);
}
double f3SecondDer(double x) {
    return 2;
}

double f5(double x) {
    if (abs(x) <= 0.333333333333333333333333) {
        return 1;
    }
    return 0;
}


double f6(double x) {
    if (abs(x) <= 0.5) {
        return 1 - 2 * abs(x);
    }
    return 0;
}

double f3Velocity(double x) {
    return cos(x);
}


double RightBoundaryEx3(double t) {
    return 2 * (t + 1);
}

double leftBoundaryEx4(double t)
{
    return sin(t);
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
    ofstream filesol(filename), filegridx("gridx.txt"), filegridt("gridt.txt");


    // делаем сетку
    vector<double> gridt{0}, gridx{0};
    filegridt << 0 << " ";
    filegridx << 0 << " ";
    vector<double> prev{},cur{}, next{};
    double t = 0, x=0;
        
    while (abs(t - T) >= 1e-10) {
        t += tau;
        gridt.emplace_back(t);
        filegridt << t << " ";
    }
    while (abs(x-L) >= 1e-10) {
        x += h;
        gridx.emplace_back(x);
        filegridx << x << " ";
    }


    int n = gridx.size(), m = gridt.size();
    // Заполняем 0ой слой и 1ый
    for (int i = 0; i < n; i++) {
        prev.push_back(InitialPositionString(gridx[i]));
        filesol << prev[i] << " ";
    }
    cur.push_back(leftBoundary(gridt[1]));
    for (int i = 1; i < n - 1; i++) {
        cur.push_back(prev[i] + InitialVelocityString(gridx[i]) * tau + a * a * tau * tau * 0.5 * InitialPositionStringSecondDerivative(gridx[i]));
    }
    cur.push_back(RightBoundary(gridt[1]));
    filesol << endl;

    if (step_t == 1) {
        for (int i = 0; i < n; i++) {
            filesol << cur[i] << " ";
        }
        filesol << endl;
    }


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
        if (j % step_t==0) {
            for (int i = 0; i < n; i++) {
                filesol << cur[i] << " ";
            }
            filesol << endl;
        }
        prev = cur;
        cur = next;
    }
}


void ExplicitMethodKolebaniaMatrix(string filename, double L, double T, double tau, double h, double a, double(&InitialPositionString)(double x),
    double(&InitialPositionStringSecondDerivative)(double x), double(&InitialVelocityString)(double x),
    double(&leftBoundary)(double x), double(&RightBoundary)(double x)) {
    ofstream filesol(filename), filegridx("gridx.txt"), filegridt("gridt.txt");
    vector<vector<double>> solMatrix{};
    vector<double> solvec{};

    // делаем сетку
    vector<double> gridt{ 0 }, gridx{ 0 };
    filegridt << 0 << " ";
    filegridx << 0 << " ";
    vector<double> prev{}, cur{}, next{};
    double t = 0, x = 0;

    while (1==1) {
        t += tau;
        if (t > T) { break; }
        gridt.emplace_back(t);
        filegridt << t << " ";
    }
    while (1==1) {
        x += h;
        if (x > L) { break; }
        gridx.emplace_back(x);
        filegridx << x << " ";
    }

    int n = gridx.size(), m = gridt.size();
    // Заполняем 0ой слой и 1ый
    for (int i = 0; i < n; i++) {
        prev.push_back(InitialPositionString(gridx[i]));
    }
    solMatrix.push_back(prev);
    cur.push_back(leftBoundary(gridt[1]));
    for (int i = 1; i < n - 1; i++) {
        cur.push_back(prev[i] + InitialVelocityString(gridx[i]) * tau + a * a * tau * tau * 0.5 * InitialPositionStringSecondDerivative(gridx[i]));
    }
    cur.push_back(RightBoundary(gridt[1]));
    solMatrix.push_back(cur);



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
        solMatrix.push_back(next);
        prev = cur;
        cur = next;
    }
    for (int i = 0; i < solMatrix.size(); i += step_t) {
        for (int j = 0; j < solMatrix[0].size(); j++) {
            filesol << solMatrix[i][j] << " ";
        }
        filesol << endl;
    }
}

// L - правая граница, left - левая граница
void ExplicitMethodKolebaniaForLeftnonzero(string filename, double L,double left, double T, double tau, double h, double a, double(&InitialPositionString)(double x),
    double(&InitialPositionStringSecondDerivative)(double x), double(&InitialVelocityString)(double x),
    double(&leftBoundary)(double x), double(&RightBoundary)(double x)) {
    ofstream filesol(filename), filegridx("gridx.txt"), filegridt("gridt.txt");
    vector<vector<double>> solMatrix{};
    vector<double> solvec{};

    // делаем сетку
    vector<double> gridt{ 0 }, gridx{left};
    filegridt << 0 << " ";
    filegridx << left << " ";
    vector<double> prev{}, cur{}, next{};
    double t = 0, x = left;

    while (1 == 1) {
        t += tau;
        if (t > T) { break; }
        gridt.emplace_back(t);
        filegridt << t << " ";
    }
    while (1 == 1) {
        x += h;
        if (x > L) { break; }
        gridx.emplace_back(x);
        filegridx << x << " ";
    }

    int n = gridx.size(), m = gridt.size();
    // Заполняем 0ой слой и 1ый
    for (int i = 0; i < n; i++) {
        prev.push_back(InitialPositionString(gridx[i]));
    }
    solMatrix.push_back(prev);
    cur.push_back(leftBoundary(gridt[1]));
    for (int i = 1; i < n - 1; i++) {
        cur.push_back(prev[i] + InitialVelocityString(gridx[i]) * tau + a * a * tau * tau * 0.5 * InitialPositionStringSecondDerivative(gridx[i]));
    }
    cur.push_back(RightBoundary(gridt[1]));
    solMatrix.push_back(cur);



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
        solMatrix.push_back(next);
        prev = cur;
        cur = next;
    }
    for (int i = 0; i < solMatrix.size(); i += step_t) {
        for (int j = 0; j < solMatrix[0].size(); j++) {
            filesol << solMatrix[i][j] << " ";
        }
        filesol << endl;
    }
}

void ExplicitMethodKolebaniaForLeftnonzeroNumDer(string filename, double L, double left, double T, double tau, double h, double a, double(&InitialPositionString)(double x), double(&InitialVelocityString)(double x),
    double(&leftBoundary)(double x), double(&RightBoundary)(double x)) {
    ofstream filesol(filename), filegridx("gridx.txt"), filegridt("gridt.txt");
    vector<vector<double>> solMatrix{};
    vector<double> solvec{};

    // делаем сетку
    vector<double> gridt{ 0 }, gridx{ left };
    filegridt << 0 << " ";
    filegridx << left << " ";
    vector<double> prev{}, cur{}, next{};
    double t = 0, x = left;

    while (t < T) {
        t += tau;
        gridt.emplace_back(t);
        filegridt << t << " ";

    }
    while (x < L) {
        x += h;
        gridx.emplace_back(x);
        filegridx << x << " ";

    }

    //DisplayVector(gridx);
    int n = gridx.size(), m = gridt.size();
    // Заполняем 0ой слой и 1ый
    for (int i = 0; i < n; i++) {
        prev.push_back(InitialPositionString(gridx[i]));
    }
    solMatrix.push_back(prev);
    cur.push_back(leftBoundary(gridt[1]));
    for (int i = 1; i < n - 1; i++) {
        cur.push_back(prev[i] + InitialVelocityString(gridx[i]) * tau + a * a * tau * tau * 0.5 * SecondDerNum(InitialPositionString,gridx[i]));
        cout << SecondDerNum(InitialPositionString, gridx[i]) << " ";
    }
    cout << endl;
    cur.push_back(RightBoundary(gridt[1]));
    DisplayVector(gridx);
    DisplayVector(prev);
    DisplayVector(cur);
    solMatrix.push_back(cur);



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
        solMatrix.push_back(next);
        prev = cur;
        cur = next;
        DisplayVector(next);
    }
    for (int i = 0; i < solMatrix.size(); i += step_t) {
        for (int j = 0; j < solMatrix[0].size(); j++) {
            filesol << solMatrix[i][j] << " ";
        }
        filesol << endl;
    }
}



void GetSolExample2(double L, double T, double tau, double h, double epsilon ) {
    ofstream fileExactSol("ExactSolExample2.txt");

    vector<double> gridt{ 0 }, gridx{ 0 };
    double t = 0, x = 0;
    while (abs(t - T) >= 1e-10) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (abs(x - L) >= 1e-10) {
        x += h;
        gridx.emplace_back(x);
    }
    int kmax = int(1 / 2 * (sqrt(2 / (pi * pi * epsilon)) - 1)) + 1;
    double curSol = 0;
    for (int i = 0; i < gridt.size(); i++) {

        for (int j = 0; j < gridx.size(); j++) {
            curSol = 0;
            for (int k = 0; k < kmax; k++) {

                curSol += 1 / ((2 * k + 1) * (2 * k + 1) * (2 * k + 1)) * sin((2 * k + 1) * pi * gridx[j]) * cos((2 * k + 1) * pi * gridt[i]);
            }
            curSol *= 8 / (pi * pi * pi);
            fileExactSol << curSol << " ";
        }
        fileExactSol << endl;
    }
}


double inftyNorm(vector<double> vec1, vector<double> vec2) {
    double res = 0;
    int n = vec1.size();
    for (int i = 0; i < n; i++) {
        res = max(res, abs(vec1[i] - vec2[i]));
    }
    return res;
}

vector<double> read_file(string file, int time_zone, int n) {
    ifstream f;
    vector<double> res;
    string temp;
    f.open(file);
    int counter = 0;
    while (counter < time_zone) {
        /*getline(f, temp);*/
        for (int i = 0; i < n + 1; i++) {
            double cur_T;
            f >> cur_T;
            //res.push_back(cur_T);
        }
        counter++;
    }
    for (int i = 0; i < n + 1; i++) {
        double cur_T;
        f >> cur_T;

        res.push_back(cur_T);
    }
    f.close();
    return res;
}

void GetOrder_p(double L, double T, double tau, double h, double a, double(&InitialPositionString)(double x),
    double(&InitialPositionStringSecondDerivative)(double x), double(&InitialVelocityString)(double x),
    double(&leftBoundary)(double x), double(&RightBoundary)(double x), double(*analytic_sol)(double, double)) {
    double sigma = 1;
    //cout <<"ust" << 0.5 - (c * rho * h * h) / (4 * tau * 1) << endl;

    ExplicitMethodKolebania("order1_h1.txt",L, T, tau, h, a, InitialPositionString, InitialPositionStringSecondDerivative,
        InitialVelocityString, leftBoundary, RightBoundary );
    ExplicitMethodKolebania("order1_h2.txt", L, T, tau, h/2, a, InitialPositionString, InitialPositionStringSecondDerivative,
        InitialVelocityString, leftBoundary, RightBoundary);
    ExplicitMethodKolebania("order1_h3.txt", L, T, tau, h/4, a, InitialPositionString, InitialPositionStringSecondDerivative,
        InitialVelocityString, leftBoundary, RightBoundary);
    vector<double> sol1, sol2, sol3;
    vector<double> gridh = { }, gridh2 = { }, gridh4 = { };
    //int tz_dynamic = 1 / tau;
    sol1 = read_file("order1_h1.txt", tz, L / h);
    /*DisplayVector(sol1);
    cout << endl;*/
    sol2 = read_file("order1_h2.txt", tz, L / (0.5 * h));
    /*DisplayVector(sol2);
    cout << endl;*/

    sol3 = read_file("order1_h3.txt", tz, L / (0.25 * h));
    /*cout << endl;
    DisplayVector(sol3);*/

    vector<double> asol1 = { }, asol2 = { }, asol3 = {  };

    double x = 0;

    while (L - x > 1e-15) {
        asol1.emplace_back(analytic_sol(x, tau * tz));
        gridh.emplace_back(x);
        x += h;
    }
    x = 0;
    while (L - x > 1e-15) {
        asol2.emplace_back(analytic_sol(x, tau * tz));
        gridh2.emplace_back(x);
        x += h / 2;
    }
    x = 0;
    while (L - x > 1e-15) {

        asol3.emplace_back(analytic_sol(x, tau * tz));
        gridh4.emplace_back(x);
        x += h / 4;
    }


    double Eh4, Eh, Eh2;
    Eh = inftyNorm(asol1, sol1);
    Eh2 = inftyNorm(asol2, sol2);
    Eh4 = inftyNorm(asol3, sol3);
    cout << setprecision(16) << "Eh1 = " << Eh << endl;
    cout << setprecision(16) << "Eh2 = " << Eh2 << endl;
    cout << setprecision(16) << "Eh4 = " << Eh4 << endl;
    double R = (Eh4 - Eh) / (Eh4 - Eh2);
    cout << log2(R - 1);

}



void GetOrder_tau(double L, double T, double tau, double h, double a, double(&InitialPositionString)(double x),
    double(&InitialPositionStringSecondDerivative)(double x), double(&InitialVelocityString)(double x),
    double(&leftBoundary)(double x), double(&RightBoundary)(double x), double(*analytic_sol)(double, double)) {
    double sigma = 1;
    //cout <<"ust" << 0.5 - (c * rho * h * h) / (4 * tau * 1) << endl;

    ExplicitMethodKolebania("order1_tau1.txt", L, T, tau, h, a, InitialPositionString, InitialPositionStringSecondDerivative,
        InitialVelocityString, leftBoundary, RightBoundary);
    ExplicitMethodKolebania("order1_tau2.txt", L, T, tau/2, h, a, InitialPositionString, InitialPositionStringSecondDerivative,
        InitialVelocityString, leftBoundary, RightBoundary);
    ExplicitMethodKolebania("order1_tau3.txt", L, T, tau/4, h, a, InitialPositionString, InitialPositionStringSecondDerivative,
        InitialVelocityString, leftBoundary, RightBoundary);
    vector<double> sol1, sol2, sol3;
    vector<double> gridh = { }, gridh2 = { }, gridh4 = { };

    sol1 = read_file("order1_tau1.txt", tz, L / h);
    /*DisplayVector(sol1);
    cout << endl;*/
    sol2 = read_file("order1_tau2.txt", tz*2, L / h);
    /*DisplayVector(sol2);
    cout << endl;*/

    sol3 = read_file("order1_tau3.txt", tz*4, L / h);
    /*cout << endl;
    DisplayVector(sol3);*/

    vector<double> asol1 = { }, asol2 = { }, asol3 = {  };

    double x = 0;

    while (L - x > 1e-15) {
        asol1.emplace_back(analytic_sol(x, tau * tz));
        gridh.emplace_back(x);
        x += h;
    }
    x = 0;

    double Etau4, Etau, Etau2;
    Etau = inftyNorm(asol1, sol1);
    Etau2 = inftyNorm(asol1, sol2);
    Etau4 = inftyNorm(asol1, sol3);
    cout << setprecision(16) << "Etau = " << Etau << endl;
    cout << setprecision(16) << "Etau2 = " << Etau2 << endl;
    cout << setprecision(16) << "Etau4 = " << Etau4 << endl;
    double R = (Etau4 - Etau) / (Etau4 - Etau2);
    cout << log2(R - 1);

}



int main()
{
    // Пример 1
    //double T = 10, L = 5, tau = 0.01, h = 0.01, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebania("Example1.txt", L, T, tau, h, a, f1, f1SecondDer, g1, g1, g1);


    // Пример 2
    //double T = 3, L = 1, tau = 0.01, h = 0.01, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebania("Example2.txt", L, T, tau, h, a, f2, f2SecondDer, g1, g1, g1);
    //GetSolExample2(L, T, tau, h, epsilon);

    // Варииант Пример 3
    //double T = 100, L = 1, tau = 0.001, h = 0.1, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebania("Example3.txt", L, T, tau, h, a, f3, f3SecondDer, f3Velocity, g1, RightBoundaryEx3);

    //Пример4 Задача 3
    //double T = 100, L = 4*pi, tau = 0.01, h = 0.1, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebaniaMatrix("Task3_Curant0.1ExactDer.txt", L, T, tau, h, a, g1, g1, g1, leftBoundaryEx4, g1);
    //Пример4 Задача 3 Числ произв
    //double T = 100, L = 4*pi, tau = 0.01, h = 0.1, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebaniaForLeftnonzeroNumDer("Task3_Curant0.1NumDer.txt", L,0, T, tau, h, a,  g1, g1, leftBoundaryEx4, g1);


    // Пример5 Задача 1

    //double T = 10,left= -2 , L = 2, tau = 0.1, h = 0.1, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebaniaForLeftnonzero("Task1_Curant1ExactDer.txt", L, left, T, tau, h, a, f5, g1, g1,g1, g1);
    
    // Пример5 Задача 1 Числ произв
    double T = 10,left= -2 , L = 2, tau = 0.1, h = 0.1, a = 1, epsilon = 1e-5;
    ExplicitMethodKolebaniaForLeftnonzeroNumDer("Task1_Curant1NumDer.txt", L, left, T, tau, h, a, f5, g1,g1, g1);

    // Пример 6 Задача2
    //double T = 10, left = -1, L = 1, tau = 0.1, h = 0.1, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebaniaForLeftnonzero("Task2_Curant1ExactDer.txt", L, left, T, tau, h, a, g1, f6, g1, g1, g1);
    //cout << SecondDerNum(f6, 0);
    // Пример 6. Задача 2. Численная производная
    //double T = 10, left = -1, L = 1, tau = 0.01, h = 0.1, a = 1, epsilon = 1e-5;
    //ExplicitMethodKolebaniaForLeftnonzeroNumDer("Task2_Curant0.1NumDer.txt", L, left, T, tau, h, a, g1, f6, g1, g1);


    // Порядки h
    //double L = 1, tau = 0.1 / (5 *2*2*2*2), h = 0.1 / (2*2*2*2), a = 1, epsilon = 1e-5;
    //double T = 1.5;
    //double T = tau * tz + 0.1;
    //GetOrder_p(L, T, tau, h, a, f1, f1SecondDer, g1, g1, g1, analsolEx1);

    //double L = 1, tau = 0.1 / (5 ), h = 0.1 / (1), a = 1, epsilon = 1e-5, T = tau * tz + 0.1;
    //GetOrder_tau(L, T, tau, h, a, f1, f1SecondDer, g1, g1, g1, analsolEx1);



}

