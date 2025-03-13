#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <vector>

using namespace std;

const double c = 460;
const double rho = 7850;
const int step_t = 15;
const int step_x = 1;

//double K(double x) {
//    if (x < 2.5) return 250;
//    else if (x < 7.5) return 250 * (x - 7.5) * 0.2 + 500 * (x - 2.5) * 0.2;
//    return 500;
//}


double K(double x) {
    return 500;
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
    double(&p1)(double t), double(&p2)(double t), double(&K)(double y), double sigma = 0.5) {
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
                    left = tau * aleft;
                    center = -tau * (aright + aleft) - h * h * c * rho;
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
                    f.push_back(-h * h * c * rho * yprev[i]);
                    progoncoefs.push_back({ left, center, right });

                }
                // Ищем коэффициент последнего уравнения
                an = 0.5 * (K(yprevs[n - 1]) + K(yprevs[n - 2]));
                left = 0;
                center = 1;
                right = 0;
                f.push_back(p1(gridt[time])/alpha1);
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
}


void SolveTempEq(string file, double tau, double h, double T, double X, double(&initial_cond)(double x), vector<double> boundary_conditions,
    double(&p1)(double t), double(&p2)(double t), double(&K)(double x), double sigma = 0.5) {
    ofstream out(file);
    ofstream out2("gridFile.txt");
    vector<double> differenceSquare{};
    double t = 0;
    double x = 0;
    vector<double> gridt = { 0 };
    vector<double> gridx = { 0 };
    double alpha1 = boundary_conditions[0];
    double alpha2 = boundary_conditions[1];
    double beta1 = boundary_conditions[2];
    double beta2 = boundary_conditions[3];
    while (t <= T) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (x < X) {
        x += h;
        gridx.emplace_back(x);
    }

    double trueSquare = (initial_cond(0) + initial_cond(gridx.back())) / 2 * gridx.back();

    for (int i = 0; i < gridx.size(); i += step_x) {
        out2 << gridx[i] << " ";
    }

    vector<double> y_next;
    vector<double> y_prev;

    int n = gridx.size();
    int m = gridt.size();

    for (int i = 0; i < n; i++) {
        y_prev.emplace_back(initial_cond(gridx[i]));
    }

    for (int i = 0; i < n; i += step_x) {
        out << y_prev[i] << " ";
    }
    out << endl;

    double left, center, right;
    double curKright = K(0.5 * h), curKleft;

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
                curKleft = K(gridx[i] - 0.5 * h);
                curKright = K(gridx[i] + 0.5 * h);;
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
            curKleft = K(gridx[n - 1] - 0.5 * h);
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
                    out << y_next[i] << " ";
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
                curKleft = K(gridx[i] - 0.5 * h);
                curKright = K(gridx[i] + 0.5 * h);;
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
            curKleft = K(gridx[n - 1] - 0.5 * h);
            // считаем прогоночный коэф-т для правого ГУ 
            left = 0;
            center = 1;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(p2(gridt[j]) /alpha2);
            y_next = Progonka(progonka_coef, f);
            if (j % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << y_next[i] << " ";
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
                curKleft = K(gridx[i] - 0.5 * h);
                curKright = K(gridx[i] + 0.5 * h);;
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
            curKleft = K(gridx[n - 1] - 0.5 * h);
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
                    out << y_next[i] << " ";
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
                curKleft = K(gridx[i] - 0.5 * h);
                curKright = K(gridx[i] + 0.5 * h);;
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
            curKleft = K(gridx[n - 1] - 0.5 * h);
            // считаем прогоночный коэф-т для правого ГУ 
            left = 0;
            center = 1;
            right = 0;
            progonka_coef.push_back({ left, center, right });
            f.push_back(p2(gridt[j]) / alpha2);
            y_next = Progonka(progonka_coef, f);

            if (j % step_t == 0) {
                for (int i = 0; i < n; i += step_x) {
                    out << y_next[i] << " ";
                }
                out << endl;
            }
            y_prev = y_next;
        }

    }
    cout << "sdfs" << endl;
    cout << Maxelement(differenceSquare) << setprecision(16);
}
    


int main()
{

    //Example1
    //SolveTempEq("test1.txt", 0.5, 0.1, 75000, 10, init_cond_test1, { 1,1,0,0 }, u1_test1, u2_test1, K, 0.5);
    //Example2
    //SolveTempEq("test2.txt", 0.5, 0.1, 50000, 10, init_cond_test1, { 0,0,-1,1 }, u1_test2, u2_test2, K, 0.5);
    //Example3 energy
    //SolveTempEq("test_energy.txt", 0.5, 0.1, 50000, 10, initial_energy, { 0,0,-1,1 }, u1_test2, u2_test2, K, 0.5);
    SolveTempEq("test1.1.txt", 0.5, 0.1, 75000, 10, init_cond_test1, { 1,1,0,0 }, u1_test1, u2_test1, K1, 0.5);
}
