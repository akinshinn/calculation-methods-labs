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


const int step_t = 1;
const int step_x = 1;

const double c = 1;
const double rho = 1;
const double pi = 3.14159265359;
const double T0 = 0.1;



double InitialConditionPorydok(double x, double y) {
    return sin(2 * pi * x) * sin(pi * y);

}

double analyticPorydok(double x, double y, double t) {
    return exp(-5 * pi * pi * t) * sin(2 * pi * x) * sin(pi * y);
}

double zero_func(double x, double y) {
    return 0;
}

double test1_func(double x, double y) {
    return 1;
}

double test2left_func(double t, double y) {
    return 1 + y;
}

double test2right_func(double t, double y) {
    return 1 + y;
}

double test2upper_func(double t, double x) {
    return 1;
}


double test2bellow_func(double t, double x) {
    return -1;
}


double test3left_func(double t, double y) {
    return 0;
}

double test3right_func(double t, double y) {
    return 2;
}

double test3upper_func(double t, double x) {
    return 1+x*x;
}

double test3bellow_func(double t, double x) {
    return x*x;
}

double answerEx3(double x, double y) {

    return x * x + y * y;
}

double integralRightSideExample3SecondTimeStep(double tau, double h1, double h2, int i) {
    return tau * h1 * h2 / 2;

}
double integralRightSideExample3FirstTimeStep(double tau, double h1, double h2, int i) {
    return tau * h1 * h2 ;

}
double integralRightSideExample4SecondTimeStep(double tau, double h1, double h2, int i) {
    return tau * h2 / 4 * (sin(2 * ((i + 0.5) * h1)) - sin(2 * i * h1));

}
double integralRightSideExample4FirstTimeStep(double tau, double h1, double h2, int i) {
    return tau * h2 / 4 * (sin(2 * ((i + 1) * h1)) - sin(2 * i * h1));

}

double integralRightSideExample2SecondTimeStep(double tau, double h1, double h2, int i) {
    return 0;

}
double integralRightSideExample2FirstTimeStep(double tau, double h1, double h2, int i) {
    return 0;

}

double test1_init(double x, double y) {
    if (x == 0 || x == 1 || y == 0 || y == 1) return 1;
    return 0;
}
double test2_init(double x, double y) {
    if (x == 0 || x == 1 || y == 0 || y == 1) return 1 + y;
    return 0;
}
double test3_init(double x, double y) {
    if (x == 0) {
        return y;
    }
    else if (y == 0) { return x * x; }
    else if ((1 - 1e-12 <= x) and (x <= 1 + 1e-12)) {
        return 2 * x + y - 1;
    }
    else if ((1-y <= 1e-14)) { return 1 + x * x; }
    else return 0;
}

double test3RightSide(double x, double y) {
    return -4;
}

double test_seminar_left(double t, double y) {
    return sin(pi*y);
}


double test_seminar_bottom(double t, double x) {
    return sin(pi*x);
}



double var1_f(double x, double y) {
    return -2 * cos(2 * x);
}

double var1_lower(double t, double x) {
    return -1;
}

double var1_upper(double t, double x) {
    return 1;
}
double var1_left(double t, double y) {
    return y;
}

double var1_right(double t, double y) {
    return 1+y;
}

double var1_init(double x, double y) {
    return 0;
}


double initial_test_seminar(double x, double y) {
    if (abs(x) < 1e-15) {
        return test_seminar_left(0, y);
    }
    if (abs(y) < 1e-15) {
        return test_seminar_bottom(0, x);
    }
    return 0;
}

double inftyNorm(vector<double> vec1, vector<double> vec2) {
    double res = 0;
    int n = vec1.size();
    for (int i = 0; i < n; i++) {
        res = max(res, abs(vec1[i] - vec2[i]));
    }
    return res;
}


double MatrixNormC(vector<vector<double>> Matrix) {
    double maxel = -1;
    for (auto& row : Matrix) {
        for (auto& el : row) {
            maxel = max(maxel, abs(el));
        }
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


void DisplayMatrix2(vector<vector<double>> Matrix) {
    for (int i = Matrix.size() - 1; i >= 0; i--) {
        for (int j = 0; j < Matrix[0].size(); j++) {
            cout << Matrix[i][j] << " ";
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




double CalculateDer_t(vector<vector<double>> M1, vector<vector<double>> M2, double tau) {
    double res = -1;
    for (int i = 0; i < M1.size(); i++) {
        for (int j = 0; j < M1[0].size(); j++) {
            double cur_der = abs(M2[i][j] - M1[i][j]) / tau;
            res = max(res, cur_der);
        }
    }
    return res;
}


void write_matrix(vector<vector<double>> M, string file) {
    ofstream out(file);
    if (out.is_open()) {
        for (int i = 0; i < M.size(); i++) {
            for (int j = 0; j < M[0].size(); j++) {
                out << M[i][j] << " " << setprecision(16);
            }
            out << endl;
        }
    }
}


vector<vector<double>> SolvePoisson(
    double T, double L1, double L2, double tau, double h1, double h2,
    double(&u_left)(double t, double y), double(&u_right)(double t, double y),
    double(&u_upper)(double t, double x), double(&u_below)(double t, double x),
    double(&f)(double x, double y), double(&initial)(double x, double y), string file,
    double eps
) {
    vector<double> gridx = { 0 }, gridy = { 0 }, gridt = { 0 };
    double t = 0, x = 0, y = 0;
    while (T - t > 5e-15) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size();
    int m = gridy.size();

    vector<vector<double>> u_next, u_prev(m, vector<double>(n, 0)), u_prev_prev;

    for (int  j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            u_prev[j][i] = initial(gridx[i], gridy[j]);
        }
    }
    //DisplayMatrix(u_prev);
    vector<vector<double>> progonka_coef;
    vector<double> d;
    for (int k = 1; k < gridt.size()-1; k++) {
        t = gridt[k];
        double t_half = t + tau / 2;
        double t_next = gridt[k+1];
        u_next = {};
        u_next.push_back(vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            u_next[0][i] = u_below(t_half, gridx[i]);
        }
        //DisplayMatrix(u_next);


        for (int j = 1; j < m-1; j++) {
            progonka_coef = {};
            d = {};
            y = gridy[j];
            progonka_coef.push_back({0,1,0});
            d.push_back(u_left(t_half, y));
            for (int i = 1; i < n-1; i++) {
                x = gridx[i];
                progonka_coef.push_back({ -1/(h1*h1), 2 / tau + 2 / (h1*h1), -1 / (h1 * h1) });
                d.push_back(2/tau * u_prev[j][i] + 1/(h2*h2)*(u_prev[j+1][i] -2 * u_prev[j][i] + u_prev[j-1][i]) + f(x, y));

            }
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_right(t_half, y));

            u_next.push_back(Progonka(progonka_coef, d));

        }
        u_next.push_back(vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            u_next[m - 1][i] = u_upper(t_half, gridx[i]);
        }
        cout << endl;
        // Целый шаг
        u_prev_prev = u_prev;
        u_prev = u_next;

        u_next = {};
        for (int j = 0; j < m; j++) {
            u_next.push_back(vector<double>(n, 0));
        }
        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][0] = u_left(t_next, gridy[j]);
        }
        //DisplayMatrix(u_next);
        for (int i = 1; i < n - 1; i++) {
            progonka_coef = {};
            d = {};
            x = gridx[i];
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_below(t_next, x));
            for (int j = 1; j < m - 1; j++) {
                y = gridy[j];
                progonka_coef.push_back({ -1 / (h2 * h2), 2 / tau + 2 / (h2 * h2), -1 / (h2 * h2) });
                d.push_back(2/tau *u_prev[j][i] + 1/(h1*h1)*(u_prev[j][i+1] - 2*u_prev[j][i] + u_prev[j][i-1]) + f(x,y));
            }
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_upper(t_next, x));
            //DisplayVector(d);
            vector<double> result_progonka = Progonka(progonka_coef, d);

            //DisplayVector(result_progonka);
            for (int j = 0; j < m; j++) {
                u_next[j][i] = result_progonka[j];
            }
            //u_next.push_back();

        }

        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][n-1] = u_right(t_next, gridy[j]);
        }

        if (CalculateDer_t(u_next, u_prev_prev, tau) < eps) {
            cout << "Final t = " << t_next << endl;
            write_matrix(u_next, file);
            return u_next;
        }
        u_prev = u_next;
        //DisplayMatrix(u_next);
        
    }

}


vector<vector<double>> SolveTeploprovodnost(
    double T, double L1, double L2, double tau, double h1, double h2,
    double(&u_left)(double t, double y), double(&u_right)(double t, double y),
    double(&u_upper)(double t, double x), double(&u_below)(double t, double x),
    double(&f)(double x, double y), double(&initial)(double x, double y),
    double eps
) {
    vector<double> gridx = { 0 }, gridy = { 0 }, gridt = { 0 };
    double t = 0, x = 0, y = 0;
    while (T - t > 5e-15) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size();
    int m = gridy.size();

    vector<vector<double>> u_next, u_prev(m, vector<double>(n, 0)), u_prev_prev;

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            u_prev[j][i] = initial(gridx[i], gridy[j]);
        }
    }
    //DisplayMatrix(u_prev);
    vector<vector<double>> progonka_coef;
    vector<double> d;
    for (int k = 1; k < gridt.size() - 1; k++) {
        t = gridt[k];
        double t_half = t + tau / 2;
        double t_next = gridt[k + 1];
        u_next = {};
        u_next.push_back(vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            u_next[0][i] = u_below(t_half, gridx[i]);
        }
        //DisplayMatrix(u_next);


        for (int j = 1; j < m - 1; j++) {
            progonka_coef = {};
            d = {};
            y = gridy[j];
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_left(t_half, y));
            for (int i = 1; i < n - 1; i++) {
                x = gridx[i];
                progonka_coef.push_back({ -1 / (h1 * h1), 2 / tau + 2 / (h1 * h1), -1 / (h1 * h1) });
                d.push_back(2 / tau * u_prev[j][i] + 1 / (h2 * h2) * (u_prev[j + 1][i] - 2 * u_prev[j][i] + u_prev[j - 1][i]) + f(x, y));
                

            }
            //cout << endl;

            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_right(t_half, y));
            //cout << "progonka_coef t = " << t << endl;
            //DisplayMatrix(progonka_coef);
            //cout << "d t = " << t << endl;
            //DisplayVector(d);
            u_next.push_back(Progonka(progonka_coef, d));
            //cout << "display sol t = " << t << endl;
            //DisplayVector(Progonka(progonka_coef, d));

        }
        u_next.push_back(vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            u_next[m - 1][i] = u_upper(t_half, gridx[i]);
        }
        //cout << endl;
        //cout << "matrix t = " << t << endl;
        //DisplayMatrix(u_next);
        // Целый шаг
        u_prev_prev = u_prev;
        u_prev = u_next;

        u_next = {};
        for (int j = 0; j < m; j++) {
            u_next.push_back(vector<double>(n, 0));
        }
        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][0] = u_left(t_next, gridy[j]);
        }
        //DisplayMatrix(u_next);
        for (int i = 1; i < n - 1; i++) {
            progonka_coef = {};
            d = {};
            x = gridx[i];
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_below(t_next, x));
            for (int j = 1; j < m - 1; j++) {
                y = gridy[j];
                progonka_coef.push_back({ -1 / (h2 * h2), 2 / tau + 2 / (h2 * h2), -1 / (h2 * h2) });
                d.push_back(2 / tau * u_prev[j][i] + 1 / (h1 * h1) * (u_prev[j][i + 1] - 2 * u_prev[j][i] + u_prev[j][i - 1]) + f(x, y));

                //progonka_coef.push_back({ 1 / (h2 * h2), -2 / tau - 2 / (h2 * h2), 1 / (h2 * h2) });
                //d.push_back(-2 / tau * u_prev[j][i] - 1 / (h1 * h1) * (u_prev[j][i + 1] - 2 * u_prev[j][i] + u_prev[j][i - 1]) - f(x, y));
            }
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_upper(t_next, x));
            //DisplayVector(d);
            vector<double> result_progonka = Progonka(progonka_coef, d);

            //DisplayVector(result_progonka);
            for (int j = 0; j < m; j++) {
                u_next[j][i] = result_progonka[j];
            }
            //u_next.push_back();

        }

        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][n - 1] = u_right(t_next, gridy[j]);
        }

        if ( abs(t - T0) < eps) {
            cout << t << endl;
            return u_next;
        }
        u_prev = u_next;
        //DisplayMatrix(u_next);

    }

}

vector<vector<double>> CountError(vector<vector<double>> chislenoeSol, double analyticsol(double x, double y, double t),
    double t, double h1, double h2, double L1, double L2) {
    vector<double> gridx = { 0 }, gridy = { 0 };
    double x = 0, y = 0;

    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size();
    int m = gridy.size();

    vector<vector<double>> error(chislenoeSol.size(), vector<double>(chislenoeSol[0].size(), 0));
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            error[j][i] = abs(chislenoeSol[j][i] - analyticsol(gridx[i], gridy[j], t));
        }
    }
    return error;
}

void countPorydok(double T, double L1, double L2,
    double(&u_left)(double t, double y), double(&u_right)(double t, double y),
    double(&u_upper)(double t, double x), double(&u_below)(double t, double x),
    double(&f)(double x, double y), double(&initial)(double x, double y),
    double eps, double analyticsol(double x, double y, double t)) {
    vector<vector<double>> sol{}, preverror{}, nexterror{};
    double h1 = 0.1, h2=0.1, tau = 0.001;
    sol = SolveTeploprovodnost(T, L1, L2, tau, h1, h2, u_left, u_right, u_upper, u_below, f, initial, eps);
    preverror = CountError(sol, analyticsol, T0, h1, h2, L1, L2);
    //cout << "sol" << endl;
    //DisplayMatrix(sol);
    //cout << endl;
    double prevNorm = MatrixNormC(preverror);
    double curNorm;
    cout << h1 << "& " << tau << "& " << prevNorm << "&" << "-" << "& \\\\ \\hline" << endl;
    for (int i = 0; i < 3; i++) {
        h1 /= 2;
        h2 /= 2;
        tau /= 2;
        sol = SolveTeploprovodnost(T, L1, L2, tau, h1, h2, u_left, u_right, u_upper, u_below, f, initial, eps);
        nexterror = CountError(sol, analyticsol, T0, h1, h2, L1, L2);
        curNorm = MatrixNormC(nexterror);
        cout << h1 << "& " << tau << "& " << curNorm << "& " << prevNorm/curNorm <<  " \\\\ \\hline" << endl;
        prevNorm = curNorm;
        preverror = nexterror;
    }
}

vector<vector<double>> SolvePoissonUpperAndLowerV4(
    double T, double L1, double L2, double tau, double h1, double h2,
    double(&u_left)(double t, double y), double(&u_right)(double t, double y),
    double(&u_upper)(double t, double x), double(&u_below)(double t, double x),
    double(&f)(double x, double y), double(&initial)(double x, double y), string file,
    double eps
) {
    vector<double> gridx = { 0 }, gridy = { 0 }, gridt = { 0 };
    double t = 0, x = 0, y = 0;
    while (T - t > 5e-15) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size();
    int m = gridy.size();

    vector<vector<double>> u_next, u_prev(m, vector<double>(n, 0)), u_prev_prev;

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            u_prev[j][i] = initial(gridx[i], gridy[j]);
        }
    }
    //DisplayMatrix(u_prev);
    vector<vector<double>> progonka_coef{};
    vector<double> d{};
    for (int k = 1; k < gridt.size() - 1; k++) {
        t = gridt[k];
        double t_half = t + tau / 2;
        double t_next = gridt[k + 1];
        u_next = {};
        u_next.push_back(vector<double>(n, 0));
        progonka_coef = {};
        d = {};

        // Посчитаем 0ыой слой y=0
        //progonka_coef.push_back({ 0,1,0 });
        //d.push_back(u_left(t_half, 0));
        //for (int i = 1; i < n - 1; i++) {
        //    progonka_coef.push_back({ -tau * h2 / (4 * h1), h1 * h2 / 2 + tau * h2 / (2 * h1), -tau * h2 / (4 * h1) });
        //    d.push_back(tau * h1 / 2 * ((u_prev[1][i] - u_prev[0][i]) / h2 + u_below(t, gridx[i])) + h1 * h2 / 2 * u_prev[0][i] + IntegralRightSideFirstTimeStep(tau, h1, h2, i));
        //}
        //progonka_coef.push_back({ 0,1,0 });
        //d.push_back(u_right(t_half, 0));
        //u_next[0] = Progonka(progonka_coef, d);
        u_next[0] = u_prev[0];

        //DisplayMatrix(u_next);

        for (int j = 1; j < m - 1; j++) {
            progonka_coef = {};
            d = {};
            y = gridy[j];
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_left(t_half, y));
            for (int i = 1; i < n - 1; i++) {
                x = gridx[i];
                progonka_coef.push_back({ -1 / (h1 * h1), 2 / tau + 2 / (h1 * h1), -1 / (h1 * h1) });
                d.push_back(2 / tau * u_prev[j][i] + 1 / (h2 * h2) * (u_prev[j + 1][i] - 2 * u_prev[j][i] + u_prev[j - 1][i]) + f(x, y));

            }
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_right(t_half, y));

            u_next.push_back(Progonka(progonka_coef, d));

        }

        u_next.push_back(vector<double>(n, 0));
        progonka_coef = {};
        d = {};
        // Верхний слой y = Y


        //progonka_coef.push_back({ 0,1,0 });
        //d.push_back(u_left(t_half, gridy.back()));
        //for (int i = 1; i < n - 1; i++) {
        //    progonka_coef.push_back({ -tau * h2 / (4 * h1), h1 * h2 / 2 + tau * h2 / (2 * h1), -tau * h2 / (4 * h1) });
        //    d.push_back(tau * h1 / 2 * (u_upper(t, gridx[i]) - (u_prev[1][i] - u_prev[0][i]) / h2) + h1 * h2 / 2 * u_prev[m - 1][i] + IntegralRightSideFirstTimeStep(tau, h1, h2, i));
        //}
        //progonka_coef.push_back({ 0,1,0 });
        //d.push_back(u_right(t_half, gridy.back()));
        //u_next[m - 1] = Progonka(progonka_coef, d);
        u_next[m - 1] = u_prev[m - 1];

        //cout << endl;

        // Целый шаг
        u_prev_prev = u_prev;
        u_prev = u_next;

        u_next = {};
        progonka_coef = {};
        d = {};
        for (int j = 0; j < m; j++) {
            u_next.push_back(vector<double>(n, 0));
        }
        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][0] = u_left(t_next, gridy[j]);
        }

        //DisplayMatrix(u_next);
        for (int i = 1; i < n - 1; i++) {
            progonka_coef = {};
            d = {};
            x = gridx[i];
            progonka_coef.push_back({0, h1*h2/2 + tau * h1/(2*h2), -tau * h1/(2*h2)});
            d.push_back(tau*h2/(4*h1)*(u_prev[0][i+1] - 2*u_prev[0][i] + u_prev[0][i-1]) +
                tau*h1/2 * u_below(t_next, gridx[i])+h1*h2/2*u_prev[0][i]+ f(gridx[i], 0)*h1*h2*tau/4);
            for (int j = 1; j < m - 1; j++) {
                y = gridy[j];
                progonka_coef.push_back({ -1 / (h2 * h2), 2 / tau + 2 / (h2 * h2), -1 / (h2 * h2) });
                d.push_back(2 / tau * u_prev[j][i] + 1 / (h1 * h1) * (u_prev[j][i + 1] - 2 * u_prev[j][i] + u_prev[j][i - 1]) + f(x, y));
            }

            progonka_coef.push_back({-tau * h1/(2*h2), h1*h2/2 + tau * h1/(2*h2),0});
            d.push_back(tau * h2 / (4 * h1) * (u_prev[m-1][i + 1] - 2 * u_prev[m-1][i] + u_prev[m-1][i - 1])
                +u_upper(t_next, gridx[i]) * tau*h1/2 + h1*h2/2 * u_prev[m-1][i]+ f(gridx[i], gridy.back()) * h1 * h2 * tau / 4);
            //DisplayVector(d);
            vector<double> result_progonka = Progonka(progonka_coef, d);

            //DisplayVector(result_progonka);
            for (int j = 0; j < m; j++) {
                u_next[j][i] = result_progonka[j];
            }
            //u_next.push_back();

        }

        //u_next.push_back(vector<double>(n, 0));
        for (int j = 0; j < m; j++) {
            u_next[j][n - 1] = u_right(t_next, gridy[j]);
        }

        if (CalculateDer_t(u_next, u_prev_prev, tau) < eps) {
            cout << "Final t = " << t_next << endl;
            write_matrix(u_next, file);
            return u_next;
        }
        /*if (t > T - 10) {
            write_matrix(u_next, file);
            return u_next;
        }*/
        //cout << CalculateDer_t(u_next, u_prev_prev, tau) << endl;
        u_prev = u_next;
        //DisplayMatrix(u_next);
        //cout << endl;
    }

}



vector<vector<double>> SolvePoissonLeftAndRightV2(
    double T, double L1, double L2, double tau, double h1, double h2,
    double(&u_left)(double t, double y), double(&u_right)(double t, double y),
    double(&u_upper)(double t, double x), double(&u_below)(double t, double x),
    double(&f)(double x, double y), double(&initial)(double x, double y), string file, double eps
) {
    vector<double> gridx = { 0 }, gridy = { 0 }, gridt = { 0 };
    double t = 0, x = 0, y = 0;
    while (T - t > 5e-15) {
        t += tau;
        gridt.emplace_back(t);
    }
    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size();
    int m = gridy.size();

    vector<vector<double>> u_next, u_prev(m, vector<double>(n, 0)), u_prev_prev;

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            u_prev[j][i] = initial(gridx[i], gridy[j]);
        }
    }
    //DisplayMatrix(u_prev);
    vector<vector<double>> progonka_coef{};
    vector<double> d{};
    for (int k = 1; k < gridt.size() - 1; k++) {
        t = gridt[k];
        double t_half = t + tau / 2;
        double t_next = gridt[k + 1];
        u_next = {};
        u_next.push_back(vector<double>(n, 0));
        progonka_coef = {};
        d = {};
        u_next[0] = u_prev[0];
        // Посчитаем 0ыой слой y=0

        for (int i = 0; i < n; i++) {
            u_next[0][i] = u_below(t_half, gridx[i]);
        }


        //DisplayMatrix(u_next);


        for (int j = 1; j < m - 1; j++) {
            progonka_coef = {};
            d = {};
            y = gridy[j];
            //progonka_coef.push_back({0, 1/h1, -1/h1});
            progonka_coef.push_back({ 0, h1 * h2 / 2 + h2 * tau / (2 * h1), -h2 * tau / (2 * h1) });
            //d.push_back(u_left(t_half, j*h2));
            d.push_back(u_left(t_half, y) * h2 * tau / 2 + (u_prev[j + 1][0] - 2 * u_prev[j][0] + u_prev[j - 1][0]) * (h1 * tau / (4 * h2)) + f(h1/4, y)*tau*h1*h2/4 + u_prev[j][0]*h1*h2/2);

            for (int i = 1; i < n - 1; i++) {
                x = gridx[i];
                progonka_coef.push_back({ -1.0 / (h1 * h1), 2.0 / tau + 2.0 / (h1 * h1), -1.0 / (h1 * h1) });
                d.push_back(2.0 / tau * u_prev[j][i] + 1.0 / (h2 * h2) * (u_prev[j + 1][i] - 2.0 * u_prev[j][i] + u_prev[j - 1][i]) + f(x, y));

            }
            //progonka_coef.push_back({ -1/h1,1/h1, 0 });
            progonka_coef.push_back({ -tau * h2 / (2 * h1), h1 * h2 / 2 + tau * h2 / (2 * h1),0 });
            //d.push_back(u_right(t_half, j*h2));
            d.push_back((u_prev[j][n - 1] * h1 * h2 * 0.5 + u_right(t_half, y)*(tau * h2 / 2) + (u_prev[j + 1][n - 1] - 2 * u_prev[j][n - 1] + u_prev[j - 1][n - 1]) * h1 * tau * 0.25 / h2 + f(gridx.back()-h1/4, y) * h2 * tau * h1 / 4));

            u_next.push_back(Progonka(progonka_coef, d));

        }
        u_next.push_back(vector<double>(n, 0));

        u_next[m - 1] = u_prev[m - 1];
        progonka_coef = {};
        d = {};
        // Верхний слой y = Y


        for (int i = 0; i < n; i++) {
            u_next[m - 1][i] = u_upper(t_half, gridx[i]);
        }

        // Целый шаг

        u_prev_prev = u_prev;
        u_prev = u_next;

        u_next = {};
        progonka_coef = {};
        d = {};
        for (int j = 0; j < m; j++) {
            u_next.push_back(vector<double>(n, 0));
        }
        //u_next.push_back(vector<double>(n, 0));
        //u_next[0][0] = u_below(t_next, 0);
        for (int j = 1; j < m-1; j++) {
            /*u_next[j][0] = h1 * u_left(t_next, j * h2) + u_prev[j][1];*/
            u_next[j][0] = u_prev[j][0];
        }

        u_next[0][0] = u_below(t_next, 0);
        u_next[0][n - 1] = u_below(t_next, gridx[n - 1]);



        //DisplayMatrix(u_next);
        for (int i = 1; i < n - 1; i++) {
            progonka_coef = {};
            d = {};
            x = gridx[i];
            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_below(t_next, gridx[i]));
            for (int j = 1; j < m - 1; j++) {
                y = gridy[j];
                progonka_coef.push_back({ -1 / (h2 * h2), 2 / tau + 2 / (h2 * h2), -1 / (h2 * h2) });
                d.push_back(2 / tau * u_prev[j][i] + 1 / (h1 * h1) * (u_prev[j][i + 1] - 2 * u_prev[j][i] + u_prev[j][i - 1]) + f(x, y));
            }

            progonka_coef.push_back({ 0,1,0 });
            d.push_back(u_upper(t_next, gridx[i]));
            vector<double> result_progonka = Progonka(progonka_coef, d);

            //DisplayVector(result_progonka);
            for (int j = 0; j < m; j++) {
                u_next[j][i] = result_progonka[j];
            }
            //u_next.push_back();

        }

        //u_next.push_back(vector<double>(n, 0));
        

        progonka_coef = {};
        d = {};
        for (int j = 1; j < m-1; j++) {
            /*u_next[j][0] = h1 * u_left(t_next, j * h2) + u_prev[j][1];*/
            u_next[j][n-1] = u_prev[j][n-1];
        }
        u_next[m - 1][0] = u_upper(t_next, 0);
        u_next[m - 1][n - 1] = u_upper(t_next, gridx.back());

        //u_next[m - 1][n-1] = (u_upper(t_next, gridx[n-1]));

        if (CalculateDer_t(u_next, u_prev_prev, tau) < eps) {
            cout << "Final t = " << t_next << endl;
            write_matrix(u_next, file);
            return u_next;
        }
        u_prev = u_next;
        //DisplayMatrix(u_next);
        cout << endl;
    }

}



void WriteAnswer(double h1, double h2,double L1, double L2, double(&answer)(double x, double y)){
    vector<double> gridx = { 0 }, gridy = { 0 };
    double x = 0, y = 0;

    while (L1 - x >= 5e-15) {
        x += h1;
        gridx.emplace_back(x);
    }
    while (L2 - y >= 5e-15) {
        y += h2;
        gridy.emplace_back(y);
    }
    int n = gridx.size(), m = gridy.size();
    vector<vector<double>> res(m, vector<double>(n, 0));
    for (int j = 0; j < m ; j++) {
        for (int i = 0; i < n; i++) {
            res[j][i] = answer(i * h1, j * h2);
        }
    }
    cout << endl;
    cout << "Answer = " << endl;
    DisplayMatrix(res);
}



double CalculateErrorForOptimalStep(vector<vector<double>> num_matrix) {
    double error = -1;
    for (int i = 0; i < num_matrix.size(); i++) {
        for (int j = 0; j < num_matrix[0].size(); j++) {
            error = max(error, abs(num_matrix[i][j] - 1));
        }
    }
    return error;
}
int main()
{

    //TEST 1
    //double T = 100000, L1 = 1, L2 = 1, tau = 0.1, h1 = 0.1, h2 = 0.1;
    ////auto u_left = test1_func, u_right = test1_func, u_upper = test1_func, u_below = test1_func;
    //vector<vector<double>> matrix = SolvePoisson(T, L1, L2, tau, h1, h2, test1_func, test1_func, test1_func, test1_func, zero_func, test1_init, "test1.txt", 1e-4);
    //DisplayMatrix(matrix);

    //TEST2
    //double T = 1000, L1 = 1, L2 = 1, tau = 0.01, h1 = 0.1, h2 = 0.1;
    //const auto u_left = test2left_func, u_right = test2right_func, u_upper = test2upper_func, u_below = test2bellow_func;
    //vector<vector<double>> matrix = SolvePoissonUpperAndLowerV4(T, L1, L2, tau, h1, h2,
    //    test2left_func, test2right_func, test2upper_func, test2bellow_func, zero_func,
    //     zero_func, "test2.txt", 1e-8);
    //DisplayMatrix(matrix);

    //Test3
    //double T = 1000, L1 = 1, L2 = 1, tau = 0.01, h1 = 0.1, h2 = 0.1;
    //vector<vector<double>> matrix = SolvePoissonLeftAndRightV2(T, L1, L2, tau, h1, h2,
    //    test3left_func, test3right_func, test3upper_func, test3bellow_func, test3RightSide, test3_init, "test3.txt", 1e-4);
    ////DisplayMatrix(matrix);
    //WriteAnswer(h1, h2, L1, L2, answerEx3);

    // Test 4. Var 1
    //double T = 1000, L1 = pi/2, L2 = 1, tau = 0.00001, h1 = 0.1, h2 = 0.1;
    //vector<vector<double>> matrix = SolvePoissonUpperAndLowerV4(T, L1, L2, tau, h1, h2,
    //    var1_left,var1_right, var1_upper, var1_lower, var1_f, zero_func, "var1.txt", 1e-4);
    //DisplayMatrix(matrix);

    // Порядки =======================================
    //double T = 100000, L1 = 1, L2 = 1, tau = 0.001, h1 = 0.1, h2 = 0.1;
    //auto u_left = test1_func, u_right = test1_func, u_upper = test1_func, u_below = test1_func;
    //vector<vector<double>> matrix = SolvePoisson(T, L1, L2, tau, h1, h2, test1_func, test1_func, test1_func, test1_func, zero_func, test1_init, "test1_0.1_0.1.txt", 1e-4);
    //DisplayMatrix(matrix);




    // Porydok
    double T = 100000, L1 = 1, L2 = 1, eps = 1e-5;
    countPorydok(T, L1, L2, zero_func, zero_func, zero_func, zero_func, zero_func, InitialConditionPorydok, eps, analyticPorydok);


    // Сравнение погрешности
    //double T = 100000, L1 = 1, L2 = 1, tau = 0.05, h1 = 0.1, h2 = 0.1;
    //vector<vector<double>> matrix = SolvePoisson(T, L1, L2, tau, h1, h2, test1_func, test1_func, test1_func, test1_func, zero_func, test1_init, "test1.txt", 1e-4);
    //cout << "Error practical optimal step = " << CalculateErrorForOptimalStep(matrix) << endl;
    //tau = 0.041;
    //matrix = SolvePoisson(T, L1, L2, tau, h1, h2, test1_func, test1_func, test1_func, test1_func, zero_func, test1_init, "test1.txt", 1e-4);
    //cout << "Error theoretical optimal step = " << CalculateErrorForOptimalStep(matrix) << endl;
}


