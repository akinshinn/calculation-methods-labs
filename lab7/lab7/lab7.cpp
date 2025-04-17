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

const double c = 1;
const double rho = 1;

const int step_t = 1;
const int step_x = 1;

//double K(double x) {
//    if (x < 2.5) return 250;
//    else if (x < 7.5) return 250 * (x - 7.5) * 0.2 + 500 * (x - 2.5) * 0.2;
//    return 500;
//}
double analytic_sol(double x, double t) {
    return exp(-(3.14159265359) * (3.14159265359) * t) * sin((3.14159265359) * x);
}

double K_order(double x) {
    return 1;
}

double init_cond_order(double x) {
    return sin(3.14159265359 * x);
}

double inftyNorm(vector<double> vec1, vector<double> vec2) {
    double res = 0;
    int n = vec1.size();
    for (int i = 0; i < n; i++) {
        res = max(res, abs(vec1[i] - vec2[i]));
    }
    return res;
}

double K(double x) {
    return 500;
}
double Knonlin(double y) {
    return 0.5 * pow(y, 2);
}

double nonlinGu(double t) {
    //cout << sqrt(2 * 25 / 0.5) * sqrt(t) << endl;;
    return sqrt(2 * 25 / 0.5) * sqrt(t);
}
double nonlinNu(double x) {
    return 0;
}
double K1(double x) {
    if (x < 2.5)
        return 250;
    else if (x < 7.5) return 250 * (x - 7.5) / ( - 5) + 500 * (x - 2.5) / 5;
    return 500;
}

    
double u1_test1(double t) {
    return 300;
}


double u2_test1(double t) {
    return 300;
}


double init_cond_test1(double x) {
    return 300 + x * (10 - x);
}



double u1_test2(double t) {
    return 0;
}


double u2_test2(double t) {
    return 0;
}


double u2ex1(double t) {
    return 0;
}
double initial_b(double x) {
    return 20*x + 20;
}

double initial_energy(double x) {
    return 2*x + 10;
}

double CalculateSquare(double h, vector<double> y) {
    double square = 0;
    for (int i = 0; i < y.size() - 1; i++) {
        square += (y[i]+y[i+1])/2 * h;
    }
    return square;
}


double Maxelement(vector<double> vec) {
    double maxel = 0;
    for (auto& el : vec) {
        maxel = max(el, maxel);
    }
    return maxel;
}

double Minelement(vector<double> vec) {
    double maxel = 0;
    for (auto& el : vec) {
        maxel = min(el, maxel);
    }
    return maxel;
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

// boundary_conditions = {alpha1, alpha2, beta1, beta2}
void SolveTempeqNonlin(string file, double tau, double h, double T, double X, double(&initial_cond)(double x), vector<double> boundary_conditions,
    double(&p1)(double t), double(&p2)(double t), double(&K)(double y)) {
    ofstream out(file);
    ofstream out2("gridFile.txt");
    double alpha1 = boundary_conditions[0];
    double alpha2 = boundary_conditions[1];
    double beta1 = boundary_conditions[2];
    double beta2 = boundary_conditions[3];
    vector<double> gridt = { 0 };
    vector<double> gridx = { 0 };
    double t = 0;
    double x = 0;
    while (t <= T) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (x < X) {
        x += h;
        gridx.emplace_back(x);
    }
    vector<double> ynext{}, yprev{}, yprevs{}, ynexts{}, f{};
    vector<vector<double>> progoncoefs{};
    int n = gridx.size();
    int m = gridt.size();
    // 0ой слой ищем 
    for (int i = 0; i < n; i++) {
        yprev.emplace_back(initial_cond(gridx[i]));
    }
    for (int i = 0; i < n; i += step_x) {
        out << yprev[i] << " ";
    }
    out << endl;
    for (int i = 0; i < gridx.size(); i += step_x) {
        cout << "123";
        out2 << gridx[i] << " ";
    }

    if ((abs(beta1) > 1e-10) && (abs(beta2) > 1e-10)) {
        for (int time = 1; time < m; time++) {
            double left, center, right, aleft, aright, a1, an;
            yprevs = yprev;
            for (int s = 0; s < 2; s++) {
                progoncoefs = {}, f = {};
                //Ищем коэффициенты 1ого уравнения
                a1 = (K(yprevs[1]) + K(yprevs[0])) / 2;
                left = 0;
                center = c * rho * h * h / 2 - h / beta1 * alpha1 * tau + a1 * tau;
                right = -a1 * tau;
                progoncoefs.push_back({ left, center, right });
                f.push_back(c * rho * h * h / 2 * yprev[0] - h * p1(gridt[time]) * tau / beta1);
                // Ищем внутренние коэффициенты
                for (int i = 1; i < n - 1; i++) {
                    aleft = 0.5 * (K(yprevs[i]) + K(yprevs[i - 1]));
                    aright = 0.5 * (K(yprevs[i + 1]) + K(yprevs[i]));
                    left = tau * aleft;
                    center = -tau * (aright + aleft) - h * h * c * rho;
                    right = tau * aright;
                    f.push_back(-h * h * c * rho * yprev[i]);
                    progoncoefs.push_back({ left, center, right });

                }
                // Ищем коэффициент последнего уравнения
                an = 0.5 * (K(yprevs[n - 1]) + K(yprevs[n - 2]));
                left = tau * an;
                center = -an * tau - h * alpha2 * tau / beta2 - c * rho * h * h / 2;
                right = 0;
                f.push_back(-c * rho * h * h / 2 * yprev[n - 1] - h / beta2 * p2(gridt[time]) * tau);
                progoncoefs.push_back({ left, center, right });
                ynexts = Progonka(progoncoefs, f);
                yprevs = ynexts;
            }
            ynext = yprevs;

            if (time % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << ynext[i] << " ";
                }
                out << endl;
            }

            yprev = ynext;
        }

    }
    else if (abs(beta1) > 1e-10) {
        for (int time = 1; time < m; time++) {
            double left, center, right, aleft, aright, a1, an;
            yprevs = yprev;
            for (int s = 0; s < 2; s++) {
                progoncoefs = {}, f = {};
                //Ищем коэффициенты 1ого уравнения
                a1 = (K(yprevs[1]) + K(yprevs[0])) / 2;
                left = 0;
                center = c * rho * h * h / 2 - h / beta1 * alpha1 * tau + a1 * tau;
                right = -a1 * tau;
                progoncoefs.push_back({ left, center, right });
                f.push_back(c * rho * h * h / 2 * yprev[0] - h * p1(gridt[time]) * tau / beta1);
                // Ищем внутренние коэффициенты
                for (int i = 1; i < n - 1; i++) {
                    aleft = 0.5 * (K(yprevs[i]) + K(yprevs[i - 1]));
                    aright = 0.5 * (K(yprevs[i + 1]) + K(yprevs[i]));
                    //aleft = 0.5 * (K(yprevs[i] + ));
                    //aright = 0.5 * K(yprevs[i + 1]);
                    left = tau * aleft;
                    center = -tau * (aright + aleft) - h * h * c * rho;
                    right = tau * aright;
                    f.push_back(-h * h * c * rho * yprev[i]);
                    progoncoefs.push_back({ left, center, right });

                }
                // Ищем коэффициент последнего уравнения
                an = 0.5 * (K(yprevs[n - 1]) + K(yprevs[n - 2]));
                left = 0;
                center = 1;
                right = 0;
                f.push_back(p2(gridt[time]) / alpha2);
                progoncoefs.push_back({ left, center, right });
                ynexts = Progonka(progoncoefs, f);
                yprevs = ynexts;
            }
            ynext = yprevs;

            if (time % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << ynext[i] << " ";
                }
                out << endl;
            }

            yprev = ynext;
        }
    }
    else if (abs(beta2) > 1e-10) {
        for (int time = 1; time < m; time++) {
            double left, center, right, aleft, aright, a1, an;
            yprevs = yprev;
            for (int s = 0; s < 2; s++) {
                progoncoefs = {}, f = {};
                //Ищем коэффициенты 1ого уравнения
                a1 = (K(yprevs[1]) + K(yprevs[0])) / 2;
                left = 0;
                center = 1;
                right = 0;
                progoncoefs.push_back({ left, center, right });
                f.push_back(p1(gridt[time]) / alpha1);
                // Ищем внутренние коэффициенты
                for (int i = 1; i < n - 1; i++) {
                    aleft = 0.5 * (K(yprevs[i]) + K(yprevs[i - 1]));
                    aright = 0.5 * (K(yprevs[i + 1]) + K(yprevs[i]));
                    left = tau * aleft;
                    center = -tau * (aright + aleft) - h * h * c * rho;
                    right = tau * aright;
                    f.push_back(-h * h * c * rho * yprev[i]);
                    progoncoefs.push_back({ left, center, right });

                }
                // Ищем коэффициент последнего уравнения
                an = 0.5 * (K(yprevs[n - 1]) + K(yprevs[n - 2]));
                left = tau * an;
                center = -an * tau - h * alpha2 * tau / beta2 - c * rho * h * h / 2;
                right = 0;
                f.push_back(-c * rho * h * h / 2 * yprev[n - 1] - h / beta2 * p2(gridt[time]) * tau);
                progoncoefs.push_back({ left, center, right });
                ynexts = Progonka(progoncoefs, f);
                yprevs = ynexts;
            }
            ynext = yprevs;

            if (time % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << ynext[i] << " ";
                }
                out << endl;
            }

            yprev = ynext;
        }
    }

    else {
        for (int time = 1; time < m; time++) {
            double left, center, right, aleft, aright, a1, an;
            yprevs = yprev;
            for (int s = 0; s < 100; s++) {
                progoncoefs = {}, f = {};
                //Ищем коэффициенты 1ого уравнения
                a1 = (K(yprevs[1]) + K(yprevs[0])) / 2;
                left = 0;
                center = 1;
                right = 0;
                progoncoefs.push_back({ left, center, right });
                f.push_back(p1(gridt[time]) / alpha1);

                // Ищем внутренние коэффициенты
                for (int i = 1; i < n - 1; i++) {
                    aleft = 0.5 * (K(yprevs[i]) + K(yprevs[i - 1]));
                    aright = 0.5 * (K(yprevs[i + 1]) + K(yprevs[i]));
                    left = aleft / (h*h);
                    center = - (aright + aleft) / (h*h) -  c * rho / tau;
                    right = aright / (h * h);
                    f.push_back(- c * rho * yprev[i] / tau);
                    progoncoefs.push_back({ left, center, right });

                }

                // Ищем коэффициент последнего уравнения
                an = 0.5 * (K(yprevs[n - 1]) + K(yprevs[n - 2]));
                left = 0;
                center = 1;
                right = 0;
                f.push_back(p2(gridt[time])/alpha2);
                progoncoefs.push_back({ left, center, right });
                ynexts = Progonka(progoncoefs, f);
                if (inftyNorm(yprevs, ynexts) < 1e-6) break;

                yprevs = ynexts;
            }
            ynext = ynexts;

            if (time % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << ynext[i] << " ";
                }
                out << endl;
            }

            yprev = ynext;
        }
    }
}


void SolveTempEq(string file, double tau, double h, double T, double X, double(&initial_cond)(double x), vector<double> boundary_conditions,
    double(&p1)(double t), double(&p2)(double t), double(&K)(double x), double sigma = 0.5) {
    ofstream out(file);
    ofstream out2("gridx.txt");
    ofstream out3("gridt.txt");
    vector<double> differenceSquare{};
    double t = 0;
    double x = 0;
    vector<double> gridt = { 0 };
    vector<double> gridx = { 0 };
    double alpha1 = boundary_conditions[0];
    double alpha2 = boundary_conditions[1];
    double beta1 = boundary_conditions[2];
    double beta2 = boundary_conditions[3];
    while (T- t > 5e-15) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (X - x > 5e-15) {
        x += h;
        //cout << (x) << " " << setprecision(16);
        //cout << (x<X) << endl;
        gridx.emplace_back(x);
    }
    //cout << (1 < 1) << endl;
    //DisplayVector(gridx);
    double trueSquare = (initial_cond(0) + initial_cond(gridx.back())) / 2 * gridx.back();

    for (int i = 0; i < gridx.size(); i += step_x) {
        out2 << gridx[i] << " ";
    }
    for (int i = 0; i < gridt.size(); i += step_t) {
        out3 << gridt[i] << " ";
    }

    vector<double> y_next;
    vector<double> y_prev;

    int n = gridx.size();
    int m = gridt.size();

    for (int i = 0; i < n; i++) {
        y_prev.emplace_back(initial_cond(gridx[i]));
    }

    for (int i = 0; i < n; i += step_x) {
        out << setprecision(16) << y_prev[i] << " ";
    }
    out << endl;

    double left, center, right;
    double curKright = (K(gridx[0]) + K(gridx[1])) / 2, curKleft;

    vector<vector<double>> progonka_coef;
    vector<double> f;

    if ((abs(beta1) > 1e-10) && (abs(beta2) > 1e-10)) {
        for (int j = 1; j < m; j++) {
            progonka_coef = {};
            f = {};
            // считаем прогоночный коэф-т для левого ГУ 
            left = 0;
            center = -c * rho * h * h * 0.5 + tau * sigma * h * alpha1 / beta1 - tau * sigma * curKright;
            right = tau * sigma * curKright;
            progonka_coef.push_back({ left, center, right });
            f.push_back(
                -tau * (1 - sigma) * (curKright * (y_prev[1] - y_prev[0])
                    + h * (alpha1 * y_prev[0] - p1(gridt[j - 1])) / beta1)
                - 0.5 * c * rho * h * h * y_prev[0] + tau * sigma * h * p1(gridt[j]) / beta1
            );
            // считаем прогоночный коэф-т для внутр точек
            for (int i = 1; i < n - 1; i++) {
                curKleft = (K(gridx[i]) + K(gridx[i - 1])) / 2;
                curKright = (K(gridx[i]) + K(gridx[i + 1])) / 2;
                left = tau * sigma * curKleft;
                center = -c * rho * h * h - tau * sigma * (curKright + curKleft);
                right = tau * sigma * curKright;
                f.push_back(

                    -c * rho * h * h * y_prev[i]
                    - tau *
                    (
                        (1 - sigma) *
                        (
                            curKright * (y_prev[i + 1] - y_prev[i]) -
                            curKleft * (y_prev[i] - y_prev[i - 1])
                            )
                        )

                );
                progonka_coef.push_back({ left, center, right });
            }
            curKleft = (K(gridx[n - 1]) + K(gridx[n - 2])) / 2;
            // считаем прогоночный коэф-т для правого ГУ 
            left = tau * sigma * curKleft;
            center = -tau * sigma * curKleft - tau * sigma * h * alpha2 / beta2 - 0.5 * c * rho * h * h;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(
                -tau * (1 - sigma) *
                (
                    h * (p2(gridt[j - 1]) - alpha2 * y_prev[n - 1]) / beta2
                    - curKleft * (y_prev[n - 1] - y_prev[n - 2])
                    )
                - 0.5 * c * rho * h * h * y_prev[n - 1]
                - tau * sigma * h * p2(gridt[j]) / beta2
            );
            y_next = Progonka(progonka_coef, f);

            differenceSquare.push_back(abs(trueSquare - CalculateSquare(h, y_next)));

            if (j % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << setprecision(16) << y_next[i] << " ";
                }
                out << endl;
            }
            y_prev = y_next;
        }
    }
    else if (abs(beta1) > 1e-10) {
        for (int j = 1; j < m; j++) {
            progonka_coef = {};
            f = {};
            // считаем прогоночный коэф-т для левого ГУ 
            left = 0;
            center = -c * rho * h * h * 0.5 + tau * sigma * h * alpha1 / beta1 - tau * sigma * curKright;
            right = tau * sigma * curKright;
            progonka_coef.push_back({ left, center, right });
            f.push_back(
                -tau * (1 - sigma) * (curKright * (y_prev[1] - y_prev[0])
                    + h * (alpha1 * y_prev[0] - p1(gridt[j - 1])) / beta1)
                - 0.5 * c * rho * h * h * y_prev[0] + tau * sigma * h * p1(gridt[j]) / beta1
            );
            // считаем прогоночный коэф-т для внутр точек
            for (int i = 1; i < n - 1; i++) {
                curKleft = (K(gridx[i]) + K(gridx[i - 1])) / 2;
                curKright = (K(gridx[i]) + K(gridx[i + 1])) / 2;
                left = tau * sigma * curKleft;
                center = -c * rho * h * h - tau * sigma * (curKright + curKleft);
                right = tau * sigma * curKright;
                f.push_back(

                    -c * rho * h * h * y_prev[i]
                    - tau *
                    (
                        (1 - sigma) *
                        (
                            curKright * (y_prev[i + 1] - y_prev[i]) -
                            curKleft * (y_prev[i] - y_prev[i - 1])
                            )
                        )

                );
                progonka_coef.push_back({ left, center, right });
            }
            curKleft = (K(gridx[n - 1]) + K(gridx[n - 2])) / 2;
            // считаем прогоночный коэф-т для правого ГУ 
            left = 0;
            center = 1;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(p2(gridt[j]) /alpha2);
            y_next = Progonka(progonka_coef, f);
            if (j % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << setprecision(16) << y_next[i] << " ";
                }
                out << endl;
            }
            y_prev = y_next;
        }
    }
    else if (abs(beta2) > 1e-10) {
        for (int j = 1; j < m; j++) {
            progonka_coef = {};
            f = {};
            // считаем прогоночный коэф-т для левого ГУ 
            left = 0;
            center = 1;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(p1(gridt[j]) /alpha1);
            // считаем прогоночный коэф-т для внутр точек
            for (int i = 1; i < n - 1; i++) {
                curKleft = (K(gridx[i]) + K(gridx[i - 1])) / 2;
                curKright = (K(gridx[i]) + K(gridx[i + 1])) / 2;
                left = tau * sigma * curKleft;
                center = -c * rho * h * h - tau * sigma * (curKright + curKleft);
                right = tau * sigma * curKright;
                f.push_back(

                    -c * rho * h * h * y_prev[i]
                    - tau *
                    (
                        (1 - sigma) *
                        (
                            curKright * (y_prev[i + 1] - y_prev[i]) -
                            curKleft * (y_prev[i] - y_prev[i - 1])
                            )
                        )

                );
                progonka_coef.push_back({ left, center, right });
            }
            curKleft = (K(gridx[n - 1]) + K(gridx[n - 2])) / 2;
            // считаем прогоночный коэф-т для правого ГУ 
            left = tau * sigma * curKleft;
            center = -tau * sigma * curKleft - tau * sigma * h * alpha2 / beta2 - 0.5 * c * rho * h * h;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(
                -tau * (1 - sigma) *
                (
                    h * (p2(gridt[j - 1]) - alpha2 * y_prev[n - 1]) / beta2
                    - curKleft * (y_prev[n - 1] - y_prev[n - 2])
                    )
                - 0.5 * c * rho * h * h * y_prev[n - 1]
                - tau * sigma * h * p2(gridt[j]) / beta2
            );
            y_next = Progonka(progonka_coef, f);

            if (j % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << setprecision(16) << y_next[i] << " ";
                }
                out << endl;
            }
            y_prev = y_next;
        }
    }
    else {
        for (int j = 1; j < m; j++) {
            progonka_coef = {};
            f = {};
            // считаем прогоночный коэф-т для левого ГУ 
            left = 0;
            center = 1;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(p1(gridt[j]) / alpha1);
            // считаем прогоночный коэф-т для внутр точек
            for (int i = 1; i < n - 1; i++) {
                curKleft = (K(gridx[i]) + K(gridx[i - 1])) / 2;
                curKright = (K(gridx[i]) + K(gridx[i + 1])) / 2;
                left = tau * sigma * curKleft;
                center = -c * rho * h * h - tau * sigma * (curKright + curKleft);
                right = tau * sigma * curKright;
                f.push_back(

                    -c * rho * h * h * y_prev[i]
                    - tau *
                    (
                        (1 - sigma) *
                        (
                            curKright * (y_prev[i + 1] - y_prev[i]) -
                            curKleft * (y_prev[i] - y_prev[i - 1])
                        )
                    )

                );
                progonka_coef.push_back({ left, center, right });
            }
            curKleft = (K(gridx[n - 1]) + K(gridx[n - 2])) / 2;
            // считаем прогоночный коэф-т для правого ГУ 
            left = 0;
            center = 1;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(p2(gridt[j]) / alpha2);
            y_next = Progonka(progonka_coef, f);

            if (j % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << setprecision(16)<<y_next[i] << " ";
                }
                out << endl;
            }
            y_prev = y_next;
        }

    }
    //cout << Maxelement(differenceSquare) << setprecision(16);
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

void GetOrder_p(double h, double tau, double(*analytic_sol)(double, double)) {
    double sigma = 1;
    //cout <<"ust" << 0.5 - (c * rho * h * h) / (4 * tau * 1) << endl;
    SolveTempEq("order1_h1.txt", tau, h, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    SolveTempEq("order1_h2.txt", tau, h/2, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    SolveTempEq("order1_h3.txt", tau, h / 4, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    vector<double> sol1, sol2, sol3;
    vector<double> gridh = { 0 }, gridh2 = { 0 }, gridh4 = { 0 };
    int tz = 100;
    sol1 = read_file("order1_h1.txt", tz + 1, 1 / h);
    sol2 = read_file("order1_h2.txt", tz + 1, 1 / (0.5*h));
    sol3 = read_file("order1_h3.txt", tz + 1, 1 / (0.25*h));
    vector<double> asol1 = { 0 }, asol2 = { 0 }, asol3 = { 0 };

    double x = 0, X = 1;
    
    while (X-x > 1e-15) {
        x += h;
        asol1.emplace_back(analytic_sol(x, tau*tz));
        gridh.emplace_back(x);
    }
    x = 0;
    while (X - x > 1e-15) {
        x += h/2;
        asol2.emplace_back(analytic_sol(x, tau * tz));
        gridh2.emplace_back(x);
    }
    x = 0;
    while (X - x > 1e-15) {
        x += h/4;
        asol3.emplace_back(analytic_sol(x, tau * tz));
        gridh4.emplace_back(x);
    }

    double Eh4, Eh, Eh2;
    Eh = inftyNorm(asol1, sol1);
    Eh2 = inftyNorm(asol2, sol2);
    Eh4 = inftyNorm(asol3, sol3);
    cout << "Eh1 = " << Eh << endl;
    cout << "Eh2 = " << Eh2 << endl;
    cout << "Eh4 = " << Eh4 << endl;
    double R = (Eh4 - Eh) / (Eh4 - Eh2);
    cout << log2(R - 1);
}

void GetOrder_q(double h, double tau, double(*analytic_sol)(double, double)) {
    double sigma = 0.5;
    //cout <<"ust" << 0.5 - (c * rho * h * h) / (4 * tau * 1) << endl;
    SolveTempEq("order1_h1.txt", tau, h, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    SolveTempEq("order1_h2.txt", tau/2, h, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    SolveTempEq("order1_h3.txt", tau/4, h , 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    vector<double> sol1, sol2, sol3;
    vector<double> gridh = { 0 }, gridh2 = { 0 }, gridh4 = { 0 };
    int tz = 100;
    sol1 = read_file("order1_h1.txt", tz + 1, 1 / h);
    sol2 = read_file("order1_h2.txt", tz*2 + 1, 1 / h);
    sol3 = read_file("order1_h3.txt", tz*4 + 1, 1 / h);
    vector<double> asol1 = { 0 }, asol2 = { 0 }, asol3 = { 0 };

    double x = 0, X = 1;

    while (X - x > 1e-15) {
        x += h;
        asol1.emplace_back(analytic_sol(x, tau * tz));
        gridh.emplace_back(x);
    }
    x = 0;
    while (X - x > 1e-15) {
        x += h;
        asol2.emplace_back(analytic_sol(x, tau * tz));
        gridh2.emplace_back(x);
    }
    x = 0;
    while (X - x > 1e-15) {
        x += h;
        asol3.emplace_back(analytic_sol(x, tau * tz));
        gridh4.emplace_back(x);
    }

    double Eh4, Eh, Eh2;
    Eh = inftyNorm(asol1, sol1);
    Eh2 = inftyNorm(asol2, sol2);
    Eh4 = inftyNorm(asol3, sol3);
    cout << "Eh1 = " << Eh << endl;
    cout << "Eh2 = " << Eh2 << endl;
    cout << "Eh4 = " << Eh4 << endl;
    double R = (Eh4 - Eh) / (Eh4 - Eh2);
    cout << log2(R - 1);
}



double CalculateError(double h1, double tau1, double(*analytic_sol)(double, double)) {
    double sigma = 1;
    vector<double> asol1 = { 0 }, asol2 = { 0 };
    vector<double> sol1{}, sol2{};
    double tz = 100; double x = 0, X = 1;
    while (X - x > 1e-15) {
        x += h1;
        asol1.emplace_back(analytic_sol(x, tau1 * tz));
        
    }

    SolveTempEq("testeror1.txt", tau1, h1, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, sigma);
    sol1 = read_file("testeror1.txt", tz, 1 / h1);

    double eh1;
    eh1 = inftyNorm(asol1, sol1);
    return eh1;

}


void getPorydokh2tau1(double h, double tau, double(*analytic_sol)(double, double)) {

    double eh1, eh2;
    eh1 = CalculateError(h, tau, analytic_sol);
    eh2 = CalculateError(h/2, tau/4, analytic_sol);
    cout << eh1 << endl;
    cout << eh2 << endl;
    cout << log2(eh1 / eh2);
}
void getPorydokh2tau2(double h, double tau, double(*analytic_sol)(double, double)) {

    double eh1, eh2;
    eh1 = CalculateError(h, tau, analytic_sol);
    eh2 = CalculateError(h / 2, tau / 2, analytic_sol);
    cout << eh1 << endl;
    cout << eh2 << endl;
    cout << log2(eh1 / eh2);
}
void testError(string f,double tau, double(*analytic_sol)(double, double)) {
    ofstream file;

    file.open(f);
    double h = 0.5;
    while (h > 1e-6) {

        file << h << " " << setprecision(16) << CalculateError(h, tau, analytic_sol) << endl;
        h /= 2;
    }

    file.close();
}








template <typename T>
std::vector<std::vector<T>> heat_eq_spatial(T rho, T c, T k1, T k2, T x1, T x2, T Q, T t0, T L, T u0, function<T(std::vector<T>)> K, std::function<T(std::vector<T>)> initial_cond, std::function<T(std::vector<T>)> left_cond, bool stream_left, std::function<T(std::vector<T>)> right_cond, bool stream_right, T tau, T h, T sigma)
{
    size_t num_points_in_space = std::round(L / h) + 1;
    size_t num_points_over_time = std::round(t0 / tau) + 1;
    std::vector<std::vector<T>> solution(1, std::vector<T>(num_points_in_space));
    std::vector<T> A(num_points_in_space - 3);
    std::vector<T> a(num_points_in_space);
    std::vector<T> B(num_points_in_space - 3);
    std::vector<T> C(num_points_in_space - 2);
    std::vector<T> F(num_points_in_space - 2);
    solution[0][0] = initial_cond({ 0 , L, u0 });;
    for (size_t i = 1; i < num_points_in_space; i++) {
        solution[0][i] = initial_cond({ i * h , L, u0 });
        // a[i] = 2.0 * K({(i - 1)*h, L, x1, x2, k1, k2, L}) * K({i*h,L, x1, x2, k1, k2}) / (K({(i - 1)*h, L, x1, x2, k1, k2, L}) + K({i*h, L, x1, x2, k1, k2}));
        a[i] = std::sqrt(K({ (i - 1) * h, L, x1, x2, k1, k2, L }) * K({ i * h, L, x1, x2, k1, k2, L }));
        // a[i] = std::sqrt(K((i - 1)*h, x1, x2, k1, k2, L) * K(i*h, x1, x2, k1, k2, L));
    }
    for (size_t j = 1; j < num_points_over_time; j++)
    {
        for (size_t i = 1; i <= num_points_in_space - 3; i++)
        {
            A[i - 1] = (sigma / h) * a[i + 1];
            B[i - 1] = (sigma / h) * a[i + 1];
        }
        for (size_t i = 1; i <= num_points_in_space - 2; i++)
        {
            C[i - 1] = -(sigma / h * a[i] + sigma / h * a[i + 1] + c * rho * h / tau);
            F[i - 1] = -(c * rho * h / tau * solution[j - 1][i] + (1 - sigma) * (a[i + 1] * (solution[j - 1][i + 1] - solution[j - 1][i]) / h - a[i] * (solution[j - 1][i] - solution[j - 1][i - 1]) / h));
        }
        std::vector<T> internal_points = TridiagonalMatrixAlgorithm(A, C, B, F);
        std::vector<T> temp{};
        // if (stream_left)
        // {
        // 	double P_tj1 = left_cond({tau*j, Q, t0});
        // 	double P_tj = left_cond({tau*(j-1), Q, t0});
        // 	double chi = (sigma * a[1] / h) / (c * rho * h / (2 * tau) + sigma * a[1] / h);
        // 	double numerator = c * rho * solution[j-1][0] * h / (2 * tau) + sigma * P_tj1 + (1 - sigma) * (P_tj - (solution[j-1][1] - solution[j-1][0]) / h);
        // 	double denominator = c * rho * h / (2 * tau) + sigma * a[1] / h;
        // 	double mu = numerator / denominator;
        // 	temp.push_back(chi*internal_points[0]+mu);
        // }
        if (stream_left)
        {
            T P_tj1 = left_cond({ tau * j, Q, t0 });
            T P_tj = left_cond({ tau * (j - 1), Q, t0 });
            T left = ((h * h + 2 * (-1 + sigma) * tau * a[1]) * solution[j - 1][0] + 2 * tau * (h * (P_tj - P_tj * sigma + P_tj1 * sigma) - (-1 + sigma) * a[1] * solution[j - 1][1] * sigma * a[1] * internal_points[0])) / (h * h + 2 * sigma * tau * a[1]);
            temp.push_back(left);
        }

        // solution[j][0] = left_cond(j*tau);
        else
        {
            temp.push_back(left_cond({ j * tau }));
        }
        temp.insert(temp.end(), internal_points.begin(), internal_points.end());
        if (stream_right)
        {
            T P_tj1 = right_cond({ tau * j, Q, t0 });
            T P_tj = right_cond({ tau * (j - 1), Q, t0 });
            T chi = (sigma * a[num_points_in_space - 1] / h) / (c * rho * h / (2 * tau) + sigma * a[num_points_in_space - 1] / h);
            T numerator = c * rho * solution[j - 1][num_points_in_space - 1] * h / (2 * tau) + sigma * P_tj1 + (1 - sigma) * (P_tj - (solution[j - 1][num_points_in_space - 1] - solution[j - 1][num_points_in_space - 2]) / h);
            T denominator = c * rho * h / (2 * tau) + sigma * a[num_points_in_space - 1] / h;
            T mu = numerator / denominator;
            temp.push_back(chi * temp[num_points_in_space - 2] + mu);
        }
        else {
            temp.push_back(right_cond({ j * tau }));
        }
        solution.emplace_back(std::move(temp));
    }
    return solution;
}









int main()
{

    //Example1
    //SolveTempEq("test1.txt", 0.5, 0.1, 75000, 10, init_cond_test1, { 1,1,0,0 }, u1_test1, u2_test1, K, 0.5);
    //Example2
    //SolveTempEq("test2.txt", 0.5, 0.1, 50000, 10, init_cond_test1, { 0,0,-1,1 }, u1_test2, u2_test2, K, 0.5);
    //Example3 energy
    //SolveTempEq("test_energy.txt", 0.5, 0.1, 50000, 10, initial_energy, { 0,0,-1,1 }, u1_test2, u2_test2, K, 0.5);
    //SolveTempEq("test1.1.txt", 0.5, 0.1, 75000, 10, init_cond_test1, { 1,1,0,0 }, u1_test1, u2_test1, K1, 0.5);
    SolveTempeqNonlin("testnonlin.txt", 2e-4, 0.2, 1, 10, nonlinNu, { 1,1,0,0 }, nonlinGu, nonlinNu, Knonlin);

    //=====================================================
    // Порядки
    //GetOrder_p(1e-1, 1e-5, analytic_sol);
    //testError("errors.txt", 0.01, analytic_sol);
    //getPorydokh2tau1(0.1, 1e-2, analytic_sol);

    //getPorydokh2tau2(0.01, 0.1, analytic_sol);
    //SolveTempEq("h0.1,tau0.1.txt", 0.1, 0.1, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 1);
    //SolveTempEq("h0.05,tau0.025.txt", 0.025, 0.05, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 1);
    //SolveTempEq("h0.025,tau0.00625.txt", 0.00625, 0.025, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 1);
    //SolveTempEq("h0.0125,tau0.0015625.txt", 0.0015625, 0.0125, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 1);
    //SolveTempEq("h0.00625,tau0.000390625.txt", 0.000390625, 0.00625, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 1);

    //SolveTempEq("h0.1,tau0.1sigma0.5.txt", 0.1, 0.1, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0.5);
    //SolveTempEq("h0.05,tau0.05sigma0.5.txt", 0.05, 0.05, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0.5);
    //SolveTempEq("h0.025,tau0.025sigma0.5.txt", 0.025, 0.025, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0.5);
    //SolveTempEq("h0.0125,tau0.0125sigma0.5.txt", 0.0125, 0.0125, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0.5);
    //SolveTempEq("h0.00625,tau0.00625sigma0.5.txt", 0.00625, 0.00625, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0.5);

    //SolveTempEq("h0.1,tau0.001sigma0.txt", 0.001, 0.1, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0);
    //SolveTempEq("h0.05,tau0.00025sigma0.txt", 0.00025, 0.05, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0);
    //SolveTempEq("h0.025,tau0.0000625sigma0.txt", 0.0000625, 0.025, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0);
    //SolveTempEq("h0.0125,tau0,000015625sigma0.txt", 0.000015625, 0.0125, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0);
    //SolveTempEq("h0.00625,tau0.00000390625sigma0.txt", 0.00000390625, 0.00625, 1, 1, init_cond_order, { 1,1,0,0 }, u1_test2, u1_test2, K_order, 0);
}
