#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <vector>
#include <functional>
#include "Header.h"
using namespace std;

const int MAXITER = 2;
const double eps = 1e-6;
const double pi = 3.141592;

void DisplayVec(vector<double> vec) {
	for (auto& el : vec) {
		cout << fixed << setprecision(10) << el << " ";
	}
	cout << endl;
}
double inftyNorm(vector<double> vec1, vector<double> vec2) {
	double res = 0;
	int n = vec1.size();
	for (int i = 0; i < n; i++) {
		res = max(res, abs(vec1[i] - vec2[i]));
	}
	return res;
}

double K(double x, double s) {
    return 0.5 * (1 - x * cos(x * s));
}

double f_test1(double x) {
    return 0.5 * (1 + sin(x));
}
double f_test2(double x) {
    return x * x + sqrt(x);
}

double Kex1(double x, double s) {
	return exp(s + x);
}

double KexIntegral1(double x) {
	return (-1 + exp(1)) * exp(x);
}

double fex1(double x) {
	return x;
}
double analSolex1(double x) {
	return x - 2.0 / (exp(2) + 1) * exp(x);
}

double Kex2(double x, double s) {
	return x* s;
}

double fex2(double x) {
	return 2.0 / 3.0 * x;
}

double solex2(double x) {
	return x;
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


vector<double> SimpleIterationMethod(double(&K)(double x, double s), double(&f)(double x), double lambda, double a,
    double b, int N, int iter = 10) {
	vector<double> gridx = GenerateUniformGrid(a, b, N);
	//DisplayVec(gridx);
	int n = gridx.size();
	double h = gridx[1] - gridx[0];
	vector<double> u_prev(n), u_next(n);
	
	for (int i = 0; i < n; i++) {
		u_prev[i] = f(gridx[i]);
	}
	//u_next = u_prev;
	//DisplayVec(u_next);
	for (int k = 0; k < iter; k++) {
		for (int i = 0; i < n; i++) {
			double integral = 0;
			//integral += (K(gridx[i], gridx[ 1]) * u_prev[1] + K(gridx[i], gridx[0]) * u_prev[0]) * h * 0.25;
			for (int j = 0; j < n - 1; j++) {
				integral += (K(gridx[i], gridx[j + 1]) * u_prev[j + 1] + K(gridx[i], gridx[j]) * u_prev[j]) * h * 0.5;
				//integral += (K(gridx[i], gridx[j]) * u_prev[j]) * h;
			}
			//integral += (K(gridx[i], gridx[n-1]) * u_prev[n-1] + K(gridx[i], gridx[n-2]) * u_prev[n-2]) * h * 0.25;
			//cout << "integral = " << integral << " x = " <<gridx[i] <<  endl;;
			u_next[i] = f(gridx[i]) + lambda * integral;

		}
		u_prev = u_next;

		//DisplayVec(u_next);


		//if (inftyNorm(u_next, u_prev) < 1e-5) break;

	}
	return u_next;

}

double f_singular(vector<double> r) {
	double phi = atan(r[1] / r[0]);
	return sin(phi);
}


double f2_singular(vector<double> r) {
	double phi = atan(r[1] / r[0]);
	return sin(2*phi);
}


vector<double> Q_singular(vector<double> r, vector<double> rho) {
	vector<double> res(2);
	double x = r[0], y = r[1], xi = rho[0], eta = rho[1];
    res[0] = 1.0 / (2.0 * pi) * (eta - y) / (pow(x - xi, 2) + pow(y - eta, 2));
    res[1] = 1.0 / (2.0 * pi) * (x - xi) / (pow(x - xi, 2) + pow(y - eta, 2));
	return res;
}

double CalculateError2(vector<double> gridx, vector<double> u_num, double(&sol)(double x)) {
	double maxerror = -1;
	for (int i = 0; i < gridx.size(); i++) {
		maxerror = max(maxerror, abs(u_num[i] - sol(gridx[i])));
	}
	return maxerror;
}

int get_main_element_row(const vector<vector<double>>& matrix, int var_row) {
	double max = abs(matrix[var_row][var_row]);
	int max_row = var_row;
	int n = matrix.size();
	for (int i = var_row + 1; i < n; ++i) {
		if (abs(matrix[i][var_row]) > max) {
			max = abs(matrix[i][var_row]);
			max_row = i;
		}
	}
	return max_row;
}


void permutate_rows(vector<vector<double>>& matrix) {
	int n = matrix.size();
	int main_elem_row;
	vector<double> temp_vec = matrix[0];
	for (int i = 0; i < n; ++i) {

		main_elem_row = get_main_element_row(matrix, i);
		temp_vec = matrix[main_elem_row];
		matrix[main_elem_row] = matrix[i];
		matrix[i] = temp_vec;
	}
}

void triangularize(vector<vector<double>>& matrix) {
	int n = matrix.size();

	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			double temp = matrix[i][j];
			if (matrix[j][j] == 0) {
				continue;
			}
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

void DisplayMatrix(vector<vector<double>> Matrix) {
	for (auto& row : Matrix) {
		for (auto el : row) {
			cout << el << " ";
		}
		cout << endl;
	}
}
vector<double> gaussMethod(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<double> x(n);

    // Прямой ход
    for (int i = 0; i < n; i++) {
        // Поиск максимального элемента для выбора главного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Перестановка строк
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        // Приведение к треугольному виду
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}
vector<double> Gauss_method(vector<vector<double>>& matrix) {
	//cout << "before triangulirized" << endl;
	//DisplayMatrix(matrix);
	//permutate_rows(matrix);
	//cout << "after permut" << endl;
	//DisplayMatrix(matrix);
	triangularize(matrix);
	//cout << "triangulirized" << endl;
	//DisplayMatrix(matrix);
	int n = matrix.size();
	vector<double> x(n);
	if (!check_matrix(matrix)) {
		cout << "det A = 0" << endl;
		//DisplayMatrix(matrix);
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
vector<double> Singular(double(&f)(vector<double> r), int N) {
	// N - число точек
	vector<vector<double>> k(N), c(N), n(N);
	for (int i = 1; i <= N; i++) {
		k[i-1] = { cos(2 * pi* (i - 0.5)/N), sin(2 * pi* (i - 0.5)/N) };
		c[i-1] = { cos(2 * pi * (i - 1)/N), sin(2 * pi  * (i - 1)/N) };
		n[i-1] = { cos(2 * pi * (i - 0.5)/N), sin(2 * pi* (i - 0.5)/N) };
	}
	double delta_l = 2 * pi / N;
	vector<double> cur_line(N + 1);
	//DisplayMatrix(k);
	write_matrix(k, "test_k.txt");
	write_matrix(c, "test_c.txt");
	vector<vector<double>> A;
	vector<double> b(N+1);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cur_line[j] = delta_l*(Q_singular(k[i], c[j])[0] * n[i][0] + Q_singular(k[i], c[j])[1] * n[i][1]);
            //DisplayVec(Q_singular(k[i], c[j]));
            //DisplayVec(k[i]);
            //DisplayVec(c[j]);
		}
		cur_line[N] = 1;
		b[i] = f(k[i]);
		A.push_back(cur_line);
	}
	// заполняем для последнего уравнения 
	for (int j = 0; j < N; j++) {
		cur_line[j] = delta_l;
	}
	cur_line[N] = 0;
	b[N] = 0;
	A.push_back(cur_line);
	//DisplayMatrix(A);
    //DisplayVec(b);
	Matrix m(A);
	vector<double> g = gaussMethod(A, b);

	//DisplayVec(g);
	return g;
}

void write_file(string file, const vector<double>& vec) {
	int n = vec.size();

	ofstream out(file);
	for (int i = 0; i < n; i++) {
		out << vec[i] << endl;
	}
}


void Setka(double h, double a, double b) {
    double curx = a;
    vector<double> gridx{ a };
    while (b - curx > 0) {
        curx += h;
        if (curx <= b) {
            gridx.push_back(curx);
        }
    }
    if (b - gridx[gridx.size() - 1] > 1e-10) {
        gridx.push_back(b);
    }
    DisplayVec(gridx);
}

vector<double> MethodKvadratur(int kolvoUzlov, double a, double b, double lambda, double(&K)(double x, double s), double(&f)(double x)) {
    double curx = a;
    vector<double> gridx{};
    gridx = GenerateUniformGrid(a, b, kolvoUzlov);
    DisplayVec(gridx);
    double h = gridx[1] - gridx[0];

    int n = gridx.size();
    vector<vector<double>> gausMatrix(n, vector<double>(n, 0));
    vector<double> rightPart{};
    for (int i = 0; i < n; i++) {
        rightPart.push_back(f(gridx[i]));
    }

    // Делам i = 0

    gausMatrix[0][0] = 1 - lambda * h * K(gridx[0], gridx[0]) / 2.0;


    for (int j = 1; j < n - 1; j++) {
        gausMatrix[0][j] = -lambda * h * K(gridx[0], gridx[j]);
    }

    gausMatrix[0][n - 1] = -lambda * h * K(gridx[0], gridx[n - 1]) / 2.0;


    // Делаем i = 1 .. n-2

    for (int i = 1; i < n - 1; i++) {
        gausMatrix[i][0] = -lambda * h * K(gridx[i], gridx[0]) / 2.0;
        for (int j = 1; j < n - 1; j++) {
            if (i == j) {
                gausMatrix[i][j] = 1 - lambda * h * K(gridx[i], gridx[j]);
            }
            else {
                gausMatrix[i][j] = -lambda * h * K(gridx[i], gridx[j]);
            }
        }
        gausMatrix[i][n - 1] = -lambda * h * K(gridx[i], gridx[n - 1]) / 2.0;
    }


    // Делаем i = n-1 последний

    gausMatrix[n - 1][0] = -lambda * K(gridx[n - 1], gridx[0]) * h / 2.0;
    for (int j = 1; j < n - 1; j++) {
        gausMatrix[n - 1][j] = -lambda * K(gridx[n - 1], gridx[j]) * h;
    }
    gausMatrix[n - 1][n - 1] = 1 - lambda * K(gridx[n - 1], gridx[n - 1]) * h / 2.0;


    //vector<double> sol = gaussMethod(gausMatrix, rightPart);

    for (int i = 0; i < n; i++) {
        gausMatrix[i].push_back(rightPart[i]);
    }

    vector<double> sol = Gauss_method(gausMatrix);

    return sol;
}



// MethodEigens - для решения интегральных уравнений с вырождженными ядерами

vector<double> MethodEigensAnalitic(double lambda, vector<vector<double>> alpha, vector<double> beta) {
    //Находит Ci 
    int m = alpha.size();

    vector<vector<double>> gausMatrix(m, vector<double>(m, 0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j) {
                gausMatrix[i][j] = 1 - lambda * alpha[i][j];
            }
            else {
                gausMatrix[i][j] = -lambda * alpha[i][j];
            }
        }
    }
    vector<double> C = gaussMethod(gausMatrix, beta);

    return C;
}

void WriteAnswerEigensAnalitic(string filename, int kolvoUzlov, double a, double b,
    double lambda, vector<vector<double>> alpha, vector<double> beta, double(&f)(double x), vector<double(*)(double)> phi) {
    vector<double> C = MethodEigensAnalitic(lambda, alpha, beta);

    int m = phi.size();
    double curx = a;
    vector<double> gridx{};
    gridx = GenerateUniformGrid(a, b, kolvoUzlov);

    ofstream file(filename);

    double curRes = 0;
    for (int i = 0; i < gridx.size(); i++) {
        curRes = f(gridx[i]);
        for (int j = 0; j < m; j++) {
            curRes += lambda * C[j] * phi[j](gridx[i]);
        }
        file << curRes << " ";
    }
    file.close();
}


vector<double> MethodEigenNumerically(double lambda, int kolvoUzlov, double a, double b, double(&f)(double x), vector<double(*)(double)> phi, vector<double(*)(double)> psi) {
    int m = psi.size();

    // посчитаем alpha, beta
    int kolvoUzlovForIntegral = 10000;

    double curx = a;
    vector<double> gridForIntegral{};
    gridForIntegral = GenerateUniformGrid(a, b, kolvoUzlovForIntegral);
    double hForIntegrals = gridForIntegral[1] - gridForIntegral[0];


    vector<vector<double>> alpha(m, vector<double>(m, 0));
    vector<double> beta(m, 0);
    double curRes = 0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            curRes = 0;
            for (int k = 0; k < gridForIntegral.size() - 1; k++) {
                curRes += hForIntegrals / 2.0 * (psi[i](gridForIntegral[k]) * phi[j](gridForIntegral[k]) + psi[i](gridForIntegral[k + 1]) * phi[j](gridForIntegral[k + 1]));
            }
            alpha[i][j] = curRes;
        }
        curRes = 0;
        for (int k = 0; k < gridForIntegral.size() - 1; k++) {
            curRes += hForIntegrals / 2.0 * (psi[i](gridForIntegral[k]) * f(gridForIntegral[k]) + psi[i](gridForIntegral[k + 1]) * f(gridForIntegral[k + 1]));
        }
        beta[i] = curRes;
    }


    vector<vector<double>> gausMatrix(m, vector<double>(m, 0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j) {
                gausMatrix[i][j] = 1 - lambda * alpha[i][j];
            }
            else {
                gausMatrix[i][j] = -lambda * alpha[i][j];
            }
        }
    }
    vector<double> C = gaussMethod(gausMatrix, beta);

    curx = a;
    vector<double> gridx{ a };
    gridx = GenerateUniformGrid(a, b, kolvoUzlov);


    int n = gridx.size();
    vector<double> res{};
    curRes = 0;
    for (int i = 0; i < n; i++) {
        curRes = f(gridx[i]);
        for (int j = 0; j < m; j++) {
            curRes += lambda * C[j] * phi[j](gridx[i]);
        }
        res.push_back(curRes);
    }

    return res;
}


void WriteAnswerEigensNumerically(string filename, double lambda, int kolvoUzlov, double a, double b, double(&f)(double x), vector<double(*)(double)> phi, vector<double(*)(double)> psi) {
    vector<double> res = MethodEigenNumerically(lambda, kolvoUzlov, a, b, f, phi, psi);

    ofstream file(filename);
    for (int i = 0; i < res.size(); i++) {
        file << res[i] << " ";
    }
    file.close();
}


template <typename T>
std::vector<T> degenerate_kernel_fredholm_solver(T lambda, T a, T b, std::function<T(T)> f, std::vector<std::function<T(T)>> epsilons, std::vector<std::function<T(T)>> psis, int grid_size) {
    const int m = epsilons.size();
    T h = (b - a) / (grid_size - 1);
    std::vector<T> weights(grid_size, h);
    weights[0] /= 2;
    weights[grid_size - 1] /= 2;
    std::vector<std::vector<T>> alpha(m, std::vector<T>(m, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            T sum = 0.0;
            for (int k = 0; k < grid_size; ++k) {
                sum += psis[i](a + k * h) * epsilons[j](a + k * h) * weights[k];
            }
            alpha[i][j] = 1.0 - lambda * sum;
        }
    }
    std::vector<T> beta(m, 0.0);
    for (int i = 0; i < m; ++i) {
        T sum = 0.0;
        for (int k = 0; k < grid_size; ++k) {
            sum += psis[i](a + k * h) * f(a + k * h) * weights[k];
        }
        beta[i] = sum;
    }
    std::vector<T> C = gaussMethod(alpha, beta);
    std::vector<T> sol(grid_size);
    for (int i = 0; i < grid_size; ++i) {
        T sum = 0.0;
        for (int j = 0; j < m; ++j) {
            sum += C[j] * epsilons[j](a + i * h);
        }
        sol[i] = f(a + i * h) + lambda * sum;
    }
    return sol;
}


double CalculateError(int kolvoUzlovminus1, double a, double b, vector<double> numSol, double(&analSol)(double x)) {
    vector<double> gridx = GenerateUniformGrid(a, b, kolvoUzlovminus1);


    double maxerror = -100;
    int n = gridx.size();
    for (int i = 0; i < n; i++) {
        //cout << numSol[i] << endl;
        //cout << abs(analSol(gridx[i]) - numSol[i]) << endl;
        maxerror = max(maxerror, abs(analSol(gridx[i]) - numSol[i]));
    }
    return maxerror;
}


// Сходиться ли порядок квадратурных формул с порядком схемы
void task1(int kolvoUzlovminus1, double a, double b, double lambda, double(&f)(double x), double(&K)(double x, double s), double(&analSol)(double x)) {

    int iter = 0;
    vector<double> res1, res2;
    double psi1, psi2;


    // для нескольких h
    while (iter < 5) {
        res1 = MethodKvadratur(kolvoUzlovminus1, a, b, lambda, K, f);
        res2 = MethodKvadratur(kolvoUzlovminus1 * 2, a, b, lambda, K, f);

        psi1 = CalculateError(kolvoUzlovminus1, a, b, res1, analSol);
        psi2 = CalculateError(kolvoUzlovminus1 * 2, a, b, res2, analSol);

        cout << kolvoUzlovminus1 << " &" << log2(psi1 / psi2) << "\\\\ \\hline" << endl;
        iter++;
        kolvoUzlovminus1 *= 2;
    }


    //Для одного h

    /*res1 = MethodKvadratur(kolvoUzlovminus1, a, b, lambda, K, f);
    res2 = MethodKvadratur(kolvoUzlovminus1 *2, a, b, lambda, K, f);

    psi1 = CalculateError(kolvoUzlovminus1, a, b, res1, analSol);
    psi2 = CalculateError(kolvoUzlovminus1*2, a, b, res2, analSol);


    cout << psi1 << endl;
    cout << psi2 << endl;
    cout << log2(psi1 / psi2) << endl;*/

}
int factorial(int n)
{
    int res = 1, i;
    for (i = 2; i <= n; i++)
        res *= i;
    return res;
}

double phiEx1Telor2(double x) {
    return 1 + x;
}

double phiEx1Telor3(double x) {
    return 1 + x + x * x / 2.0;
}

double phiEx1Telor4(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3);
}


double phiEx1Telor5(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4);
}

double phiEx1Telor6(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5);
}

double phiEx1Telor7(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5) + pow(x, 6) / factorial(6);
}


double phiEx1Telor8(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5) + pow(x, 6) / factorial(6) + pow(x, 7) / factorial(7);
}


double phiEx1Telor9(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5) +
        pow(x, 6) / factorial(6) + pow(x, 7) / factorial(7) + pow(x, 8) / factorial(8);
}


double phiEx1Telor10(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5) +
        pow(x, 6) / factorial(6) + pow(x, 7) / factorial(7) + pow(x, 8) / factorial(8) + pow(x, 9) / factorial(9);
}

double phiEx1Telor11(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5) +
        pow(x, 6) / factorial(6) + pow(x, 7) / factorial(7) + pow(x, 8) / factorial(8) + pow(x, 9) / factorial(9) + pow(x, 10) / factorial(10);
}

double phiEx1Telor12(double x) {
    return 1 + x + x * x / 2.0 + pow(x, 3) / factorial(3) + pow(x, 4) / factorial(4) + pow(x, 5) / factorial(5) +
        pow(x, 6) / factorial(6) + pow(x, 7) / factorial(7) + pow(x, 8) / factorial(8) + pow(x, 9) / factorial(9) + pow(x, 10) / factorial(10) + pow(x, 11) / factorial(11);
}

double phiEx1(double x) {
    return exp(x);
}

void task2(double lambda, double a, double b, int kolvoUzlov, vector<double(*)(double)> phi, double(&analSol)(double x), double(&f)(double x)) {
    vector<double> res;
    vector<double> gridx = GenerateUniformGrid(a, b, kolvoUzlov);
    double error = 0;
    for (int i = 0; i < phi.size(); i++) {
        res = MethodEigenNumerically(lambda, kolvoUzlov, a, b, f, { phi[i] }, { phi[i] });

        //res = degenerate_kernel_fredholm_solver<double>(lambda, a, b, f, { phi[i] }, { phi[i] }, kolvoUzlov);

        error = CalculateError(kolvoUzlov, a, b, res, analSol);
        cout << i + 2 << " &" << error << "\\\\ \\hline" << endl;
    }
}
double anal_Sol_test1_methoda(double x) {
    return 1;
}
double K_SIM_test(double x, double s) {
    return x * s;
}
double f_SIM_test(double x) {
    return 2.0 / 3.0 * x;
}
double anal_sol_SIM_test(double x) {
    return x;
}

double K_variant(double x, double s) {
    return x * sin(x * s);
}

double f1_variant(double x) {
    return cos(x);
}

double f2_variant(double x) {
    return x * x + sqrt(x);
}

int main()
{
	//vector<double> u = SimpleIterationMethod(Kex2, fex2, 1, 0, 1, 100);
	//vector<double> u = simple_iteration_fredholm_solve<double>(1, 0, 1, Kex2, fex2, 100, 10, 1e-5);
	//DisplayVec(u);
	//cout <<"error = " << CalculateError2(GenerateUniformGrid(0, 1, u.size()-1), u, solex2);
	//write_file("res1.txt",(Singular(f2_singular, 100)));
    //vector<double> res = SimpleIterationMethod(K, f_test2, 1, 0.1, 1, 200);
    //vector<double> res = MethodKvadratur(200, 0.1, 1, 1, K_variant, f2_variant);
    //vector<double> res = MethodKvadratur(200, 0.1, 1, 1, K, f_test2);
    //vector<double> res = MethodEigenNumerically(1, 200, 0, 1, f_test1, )
    //write_file("test1_1000_lim1.txt", res);

    //int N = 20;
    //auto gridx = GenerateUniformGrid(0, 1, N);
    //vector<double> iteration;
    //for (int k = 1; k <= 20; k += 1) {
    //    cout << k << " & " << CalculateError(N, 0, 1, SimpleIterationMethod(K_SIM_test, f_SIM_test, 1, 0, 1, N, k), anal_sol_SIM_test) << "\\\\ \\hline " << endl;
    //    iteration.push_back(CalculateError(N, 0, 1, SimpleIterationMethod(K_SIM_test, f_SIM_test, 1, 0, 1, N, k), anal_sol_SIM_test));
    //}
    //write_file("iteration.txt", iteration);
    vector<double> R;
    vector<double> res_Singular = Singular(f_singular, 5);

    //for (int k = 1; k <= 20; k += 1) {
    //    res_Singular = Singular(f_singular, k);
    //    //cout << res_Singular.size() << endl;
    //    double R_c = res_Singular.back();
    //    //cout << R_c << endl;
    //    R.push_back(R_c);
    //    //cout << k << " & " << R.back() << "\\\\ \\hline " << endl;
    //}
    //DisplayVec(R);
    write_file("R.txt", R);
    //write_file("iteration.txt", iteration);
    //write_file("res1.txt", res_Singular);
}

