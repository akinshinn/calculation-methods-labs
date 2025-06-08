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
const int batch_size = 1000;


double f1_variant(double t, vector<double> vars) {
	return vars[1];
}

double f2_variant(double t, vector<double> vars) {
	return 2 * vars[1] * 0.3 * (1 - vars[0] * vars[0]) + vars[0];
}

double true_x_pendulum(double t) {
	return cos(sqrt(10) * t);
}


double true_y_pendulum(double t) {
	return -sqrt(10) * sin(sqrt(10) * t);
}

vector<double> sol_pendulum(double t) {
	return { 2 * cos(sqrt(k / m) * t), -2 * sqrt(k / m) * sin(sqrt(k / m) * t) };
}

double f1_pendulum(double t, vector<double> values) {
	return values[1];
}

double f2_pendulum(double t, vector<double> values) {
	return -k / m * values[0];
}


double f1_system1_book(double t, vector<double> values) {
	return 2 * values[0] + values[1] * values[1] - 1;
}


double f2_system1_book(double t, vector<double> values) {
	return 6 * values[0] - values[1] * values[1] + 1;
}


double f1_system2_book(double t, vector<double> values) {
	return 1 - values[0] * values[0] - values[1] * values[1];
}


double f2_system2_book(double t, vector<double> values) {
	return 2 * values[0];
}


double f1_system3_book(double t, vector<double> values) {
	return 10 * (values[1] - values[0]);
}


double f2_system3_book(double t, vector<double> values) {
	return values[0] * (28 - values[2]) - values[1];
}


double f3_system3_book(double t, vector<double> values) {
	return values[0] * values[1] - 8 * values[2] / 3;
}


double f1_system1_test(double t, vector<double> values) {
	return values[1];
}


double f2_system1_test(double t, vector<double> values) {
	return values[0];
}

double f1_system2_test(double t, vector<double> values) {
	return 4 * values[0] - values[1] + exp(3 * t) * (t + sin(t));
}


double f2_system2_test(double t, vector<double> values) {
	return values[0] + 2 * values[1] + t * exp(3 * t) * cos(t);
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
	for (int fun = 0; fun < n; fun++) {
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
		for (int var = 1; var < n + 1; var++) {
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

	/*for (int i = 1; i < 4; i++) {
		vector<double> tmp{};

		for (int j = 0; j < n; j++) {
			tmp.push_back(res[i - 1][j] + tau * f[j](time[i - 1], res[i - 1]));
		}
		res.push_back(tmp);
	}*/
	// Найдем y1, y2, y3 с помощью метода Рунге-Кутты 4-го порядка
	for (int i = 1; i < 4; i++) {
		vector<double> tmp{}, k1{}, k2{}, k3{}, k4{}, arg2{}, arg3{}, arg4{};
		// Находим k
		for (int j = 0; j < n; j++) {
			k1.push_back(f[j](time[i - 1], res[i - 1]));
			arg2.push_back(res[i - 1][j] + tau / 2 * k1[j]);
		}
		for (int j = 0; j < n; j++) {
			k2.push_back(f[j](time[i - 1] + tau / 2, arg2));
			arg3.push_back(res[i - 1][j] + tau / 2 * k2[j]);
		}
		for (int j = 0; j < n; j++) {
			k3.push_back(f[j](time[i - 1] + tau / 2, arg3));
			arg4.push_back(res[i - 1][j] + tau * k3[j]);
		}
		for (int j = 0; j < n; j++) {
			k4.push_back(f[j](time[i - 1] + tau, arg4));
		}
		for (int j = 0; j < n; j++) {
			tmp.push_back(res[i - 1][j] + tau * 1 / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
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



vector<vector<double>>	(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f) {
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

    /*for (int i = 1; i < 4; i++) {
        vector<double> tmp{};

        for (int j = 0; j < n; j++) {
            tmp.push_back(res[i - 1][j] + tau * f[j](time[i - 1], res[i - 1]));
        }
        res.push_back(tmp);
    }*/

    for (int i = 1; i < 4; i++) {
        vector<double> tmp{}, k1{}, k2{}, k3{}, k4{}, arg2{}, arg3{}, arg4{};
        // Находим k
        for (int j = 0; j < n; j++) {
            k1.push_back(f[j](time[i - 1], res[i - 1]));
            arg2.push_back(res[i - 1][j] + tau / 2 * k1[j]);
        }
        for (int j = 0; j < n; j++) {
            k2.push_back(f[j](time[i - 1] + tau / 2, arg2));
            arg3.push_back(res[i - 1][j] + tau / 2 * k2[j]);
        }
        for (int j = 0; j < n; j++) {
            k3.push_back(f[j](time[i - 1] + tau / 2, arg3));
            arg4.push_back(res[i - 1][j] + tau * k3[j]);
        }
        for (int j = 0; j < n; j++) {
            k4.push_back(f[j](time[i - 1] + tau, arg4));
        }
        for (int j = 0; j < n; j++) {
            tmp.push_back(res[i - 1][j] + tau * 1 / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
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
	//int center = sol1.size() / 3.9;
	int center = sol1.size() / 2.8;
	cout << log(abs((sol3[center / (q * q)][0] - sol2[center / q][0]) / (sol2[center / q][0] - sol1[center][0]))) / log(q);
}
void PorydokAnal(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f, double q, vector<vector<double>>(&solver)(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f), vector<double>(&anal)(double t)) {
	vector<vector<double>> sol1 = solver(tau, T, y0, f);
	vector<vector<double>> sol2 = solver(tau * q, T, y0, f);

	int center = sol1.size() /2.8;
	cout << log(abs((sol2[center / (q)][0] - anal(center*tau)[0]) / (sol1[center][0] - anal(center * tau)[0]))) / log(q);
}

void ProcessAitken(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f, double q, vector<vector<double>>(&solver)(vector<double> y0,
	const vector<double>& grid,
	vector<double(*)(double, vector<double>)> syst)) {
	vector<double> time = {}, time1 = {}, time2 = {};
	double t00 = 0;
	while (t00 <= T) {
		time.push_back(t00);
		t00 += tau;
	}
	t00 = 0;
	while (t00 <= T) {
		time1.push_back(t00);
		t00 += q * tau;
	}
	t00 = 0;
	while (t00 <= T) {
		time2.push_back(t00);
		t00 += q * q * tau;
	}

	vector<vector<double>> sol1 = solver(y0, time, f);

	vector<vector<double>> sol2 = solver(y0, time1, f);
	vector<vector<double>> sol3 = solver(y0, time2, f);
	//int center = sol1.size() / 3.9;
	int center = sol1.size() / 1000;
	cout << log(abs((sol3[center / (q * q)][0] - sol2[center / q][0]) / (sol2[center / q][0] - sol1[center][0]))) / log(q);
}


void GetOrder(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f, double q, vector<vector<double>>(&solver)(double tau, double T, vector<double> y0, vector<double(*)(double, vector<double>)> f), vector<double(*)(double)> true_sol) {
		vector<vector<double>> sol1 = solver(tau, T, y0, f);
		vector<vector<double>> sol2 = solver(tau * q, T, y0, f);
		//int center = sol1.size() / 3.9;
		int center = sol1.size() / 4;
		cout << log(abs((sol2[center / q][0] - true_sol[0](T/tau * center) / (sol1[center][0] - true_sol[0](T / tau * center))))) / log(q);
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
			y_cur[j] = res[i - 1][j] + tau * syst[j](grid[i - 1], res[i - 1]);

		}
		res.emplace_back(y_cur);
	}

	return res;
}

// last_id - последний индекс, когда нужно перестать записывать (в случае, если вектор GridFunc перезаписан не полностью)
void AddFileBatch(vector<vector<double>> GridFunc, string file, int last_id = 0) {
	int num_vars = GridFunc[0].size();
	int n;
	if (last_id) {
		n = last_id;
	}
	else {
		n = GridFunc.size();
	}
	ofstream out(file, ios::app);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < num_vars; j++) {
			out << GridFunc[i][j] << " ";
		}
		out << endl;
	}
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

void BatchedEulerExplicit(
	vector<double> y0,
	const vector<double>& grid,
	vector<double(*)(double, vector<double>)> syst,
	string file) {

	int n = grid.size(); // количество узлов
	int num_vars = syst.size(); // количество переменных
	double tau = grid[1] - grid[0];
	vector<vector<double>> res(min(batch_size, n)); // res хранит вектор по каждой слою времени, 
	// в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i
	int batch_nums = n / batch_size; // кол-во пакетов


	// Первый слой пакета делаем отдельно засчет автоматического добавления Н.У.
	res[0] = y0; // добавляем начальный слой времени

	// i - номер узла, j - номер переменной
	for (int i = 1; i < min(batch_size, n); i++) {
		vector<double> y_cur(num_vars);
		for (int j = 0; j < num_vars; j++) {
			y_cur[j] = res[i - 1][j] + tau * syst[j](grid[i - 1], res[i - 1]);

		}
		res[i] = y_cur;
	}
	AddFileBatch(res, file); // Добавление в файл res
	n -= batch_size + 1;
	// Теперь цикл по каждому пакету
	vector<double> prev_point; // значение (y1, y2, ..., yn) на предыдущем узле
	for (int k = 1; k < batch_nums; k++) {
		prev_point = res[res.size() - 1]; // в начале пред точка - последняя из предыдущего пакета

		for (int i = 0; i < batch_size; i++) {
			vector<double> y_cur(num_vars);
			for (int j = 0; j < num_vars; j++) {
				y_cur[j] = prev_point[j] + tau * syst[j](grid[i - 1 + batch_size * k], prev_point);

			}
			res[i] = y_cur;
			prev_point = y_cur;
		}
		AddFileBatch(res, file);
		n -= batch_size;
	}

	// Последняя итерация по маленькому пакету (если кол-во точек не делится на размер пакета)

	prev_point = res[res.size() - 1]; // в начале пред точка - последняя из предыдущего пакета

	for (int i = 0; i < n; i++) {
		vector<double> y_cur(num_vars);
		for (int j = 0; j < num_vars; j++) {
			y_cur[j] = prev_point[j] + tau * syst[j](grid[i - 1 + batch_size * k], prev_point);

		}
		res[i] = y_cur;
		prev_point = y_cur;
	}
	AddFileBatch(res, file, n); // последняя перезаписанная точка n-1


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
			auto cur_f = syst[j];
			auto last_t = grid[i - 1];
			auto last_y = res[i - 1];
			auto func{ [tau, i, j,last_t, cur_f ,last_y](vector<double> coordinates)
				{
					double t = coordinates[0]; // syst[j] принимает на вход t, y, поэтому выделяем t
					coordinates.erase(coordinates.begin());
					double ans = ((coordinates[j] - last_y[j]) / tau);

					ans -= 0.5 * (cur_f(last_t, last_y) + cur_f(t, coordinates));
					return ans;
				}
			};
			//int iteration = 0;
			//vector<double> prevstep = y0;
			//vector<double> prevPoint;
			//prevPoint.emplace_back(t0);
			//prevPoint.insert(prevPoint.end(), prevstep.begin(), prevstep.end());
			//vector<double> curstep(n), b(n), delta{};
			//vector<vector<double>> JacobyMatrix;
			//while (iteration < Maxiter) {
			//	for (int i = 0; i < n; i++) {
			//		b[i] = -f[i](prevPoint); // правая часть системы
			//	}


			//	JacobyMatrix = ComputeJacoby(f, prevPoint);
			//	// Сводим к виду чтобы скормить Гауссу
			//	for (int i = 0; i < n; i++) {
			//		JacobyMatrix[i].push_back(b[i]);
			//	}

			//	delta = Gauss_method(JacobyMatrix);


			//	for (int i = 0; i < n; i++) {
			//		curstep[i] = prevstep[i] + delta[i];
			//		prevPoint[i + 1] = curstep[i]; // i+1, т.к. prevPoint[0] = t0
			//	}

			//	if (inftyNorm(curstep, prevstep) < 1e-6) { break; }
			//	prevstep = curstep;
			//	iteration++;
			//}

			functions.emplace_back(func);

		}
		y_cur = NewtonMethod(grid[i], res[i - 1], functions); // в качестве времени считаем текущее grid[i], а за начальное приближение считаем предыдущее решение res[i - 1]
		res.emplace_back(y_cur);
	}
	return res;

}


void BatchedSymmetricalScheme(vector<double> y0,
	const vector<double>& grid,
	vector<double(*)(double, vector<double>)> syst, string file) {
	int n = grid.size(); // количество узлов
	int num_vars = syst.size(); // количество переменных
	double tau = grid[1] - grid[0];
	vector<vector<double>> res(min(batch_size, n)); // res хранит вектор по каждой слою времени, 
	// в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i
	int batch_nums = n / batch_size; // кол-во пакетов

	res[0] = y0;
	// i - номер узла, j - номер переменной
	for (int i = 1; i < min(batch_size, n); i++) {
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
		res[i] = y_cur;
	}
	AddFileBatch(res, file); // Добавление в файл res
	n -= batch_size + 1;


	vector<double> prev_point; // значение (y1, y2, ..., yn) на предыдущем узле
	for (int k = 1; k < batch_nums; k++) {
		prev_point = res[res.size() - 1]; // в начале пред точка - последняя из предыдущего пакета

		for (int i = 0; i < batch_size; i++) {
			vector<double> y_cur;
			vector<function<double(vector<double>)>> functions; // система нелинейных уравнений
			for (int j = 0; j < num_vars; j++) {
				auto func{ [tau, i, j,prev_point,syst, grid,k](vector<double> coordinates)
					{
						double t = coordinates[0]; // syst[j] принимает на вход t, y, поэтому выделяем t
						coordinates.erase(coordinates.begin());
						double ans = ((coordinates[j] - prev_point[j]) / tau);

						ans -= 0.5 * (syst[j](grid[i - 1 + k * batch_size], prev_point) + syst[j](t, coordinates));
						return ans;
					}
				};

				functions.emplace_back(func);

			}
			y_cur = NewtonMethod(grid[i + k * batch_size], prev_point, functions); // в качестве времени считаем текущее grid[i], а за начальное приближение считаем предыдущее решение res[i - 1]
			prev_point = y_cur;
			res[i] = y_cur;
		}
		AddFileBatch(res, file);
		n -= batch_size;

	}


	prev_point = res[res.size() - 1]; // в начале пред точка - последняя из предыдущего пакета

	for (int i = 0; i < n; i++) {
		vector<double> y_cur;
		vector<function<double(vector<double>)>> functions; // система нелинейных уравнений
		for (int j = 0; j < num_vars; j++) {
			auto func{ [tau, i, j,prev_point,syst, grid, batch_nums](vector<double> coordinates)
				{
					double t = coordinates[0]; // syst[j] принимает на вход t, y, поэтому выделяем t
					coordinates.erase(coordinates.begin());
					double ans = ((coordinates[j] - prev_point[j]) / tau);

					ans -= 0.5 * (syst[j](grid[i - 1 + batch_nums * batch_size], prev_point) + syst[j](t, coordinates));
					return ans;
				}
			};

			functions.emplace_back(func);

		}
		y_cur = NewtonMethod(grid[i + k * batch_size], prev_point, functions); // в качестве времени считаем текущее grid[i], а за начальное приближение считаем предыдущее решение res[i - 1]
		prev_point = y_cur;
		res[i] = y_cur;
	}
	AddFileBatch(res, file, n); // последняя перезаписанная точка n-1


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



vector<vector<double>> Runge_Kutta(vector<double> y0,
	const vector<double>& grid,
	vector<double(*)(double, vector<double>)> syst) {
	int n = grid.size(); // количество узлов
	int num_vars = syst.size(); // количество переменных
	double tau = grid[1] - grid[0];
	vector<vector<double>> res; // res хранит вектор по каждой слою времени, 
	// в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i

	res.emplace_back(y0);
	vector<double> k1(num_vars), k2(num_vars), k3(num_vars), k4(num_vars), K_res(num_vars);
	vector<double> y_k1(num_vars), y_k2(num_vars), y_k3(num_vars);

	for (int i = 1; i < n; i++) {
		vector<double> y_cur = res[i - 1];
		y_k1 = res[i - 1];
		y_k2 = res[i - 1];
		y_k3 = res[i - 1];
		// считаем k1
		for (int j = 0; j < num_vars; j++) {
			k1[j] = syst[j](grid[i - 1], res[i - 1]);
			y_k1[j] += tau * (k1[j]) / 2;
		}

		// считаем k2 
		for (int j = 0; j < num_vars; j++) {
			k2[j] = syst[j](grid[i - 1] + tau / 2, y_k1);
			y_k2[j] += tau * (k2[j]) / 2;
		}

		// считаем k3
		for (int j = 0; j < num_vars; j++) {
			k3[j] = syst[j](grid[i - 1] + tau / 2, y_k2);
			y_k3[j] += tau * (k3[j]);
		}

		// считаем k4
		for (int j = 0; j < num_vars; j++) {
			k4[j] = syst[j](grid[i - 1] + tau, y_k3);
			K_res[j] = (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
			y_cur[j] += tau * K_res[j];

		}
		res.emplace_back(y_cur);


	}
	return res;
}



void BatchedRungeKutta(vector<double> y0,
	const vector<double>& grid,
	vector<double(*)(double, vector<double>)> syst,
	string file) {
	int n = grid.size(); // количество узлов
	int num_vars = syst.size(); // количество переменных
	double tau = grid[1] - grid[0];
	vector<vector<double>> res(min(batch_size, n)); // res хранит вектор по каждой слою времени, 
	// в каждом слое значения переменных, т.е. res[i] -> {x1, x2,..., x_num_vars} | t->t_i
	int batch_nums = n / batch_size; // кол-во пакетов

	res[0] = y0;
	vector<double> k1(num_vars), k2(num_vars), k3(num_vars), k4(num_vars), K_res(num_vars);
	vector<double> y_k1(num_vars), y_k2(num_vars), y_k3(num_vars);

	for (int i = 1; i < min(n, batch_size); i++) {
		vector<double> y_cur = res[i - 1];
		y_k1 = res[i - 1];
		y_k2 = res[i - 1];
		y_k3 = res[i - 1];
		// считаем k1
		for (int j = 0; j < num_vars; j++) {
			k1[j] = syst[j](grid[i - 1], res[i - 1]);
			y_k1[j] += tau * (k1[j]) / 2;
		}

		// считаем k2 
		for (int j = 0; j < num_vars; j++) {
			k2[j] = syst[j](grid[i - 1] + tau / 2, y_k1);
			y_k2[j] += tau * (k2[j]) / 2;
		}

		// считаем k3
		for (int j = 0; j < num_vars; j++) {
			k3[j] = syst[j](grid[i - 1] + tau / 2, y_k2);
			y_k3[j] += tau * (k3[j]);
		}

		// считаем k4
		for (int j = 0; j < num_vars; j++) {
			k4[j] = syst[j](grid[i - 1] + tau, y_k3);
			K_res[j] = (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
			y_cur[j] += tau * K_res[j];

		}
		res[i] = y_cur;
	}
	AddFileBatch(res, file);
	n -= batch_size;




	vector<double> prev_point; // значение (y1, y2, ..., yn) на предыдущем узле
	for (int k = 1; k < batch_nums; k++) {
		prev_point = res[res.size() - 1]; // в начале пред точка - последняя из предыдущего пакета
		for (int i = 0; i < batch_size; i++) {
			vector<double> y_cur = prev_point;
			y_k1 = prev_point;
			y_k2 = prev_point;
			y_k3 = prev_point;
			// считаем k1
			for (int j = 0; j < num_vars; j++) {
				k1[j] = syst[j](grid[i - 1 + k * batch_size], prev_point);
				y_k1[j] += tau * (k1[j]) / 2;
			}

			// считаем k2 
			for (int j = 0; j < num_vars; j++) {
				k2[j] = syst[j](grid[i - 1 + k * batch_size] + tau / 2, y_k1);
				y_k2[j] += tau * (k2[j]) / 2;
			}

			// считаем k3
			for (int j = 0; j < num_vars; j++) {
				k3[j] = syst[j](grid[i - 1 + k * batch_size] + tau / 2, y_k2);
				y_k3[j] += tau * (k3[j]);
			}

			// считаем k4
			for (int j = 0; j < num_vars; j++) {
				k4[j] = syst[j](grid[i - 1 + k * batch_size] + tau, y_k3);
				K_res[j] = (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
				y_cur[j] += tau * K_res[j];

			}
			res[i] = y_cur;
		}

		AddFileBatch(res, file);
		n -= batch_size;

	}


	prev_point = res[res.size() - 1]; // в начале пред точка - последняя из предыдущего пакета

	for (int i = 0; i < n; i++) {
		vector<double> y_cur = prev_point;
		y_k1 = prev_point;
		y_k2 = prev_point;
		y_k3 = prev_point;
		// считаем k1
		for (int j = 0; j < num_vars; j++) {
			k1[j] = syst[j](grid[i - 1 + k * batch_size], prev_point);
			y_k1[j] += tau * (k1[j]) / 2;
		}

		// считаем k2 
		for (int j = 0; j < num_vars; j++) {
			k2[j] = syst[j](grid[i - 1 + k * batch_size] + tau / 2, y_k1);
			y_k2[j] += tau * (k2[j]) / 2;
		}

		// считаем k3
		for (int j = 0; j < num_vars; j++) {
			k3[j] = syst[j](grid[i - 1 + k * batch_size] + tau / 2, y_k2);
			y_k3[j] += tau * (k3[j]);
		}

		// считаем k4
		for (int j = 0; j < num_vars; j++) {
			k4[j] = syst[j](grid[i - 1 + k * batch_size] + tau, y_k3);
			K_res[j] = (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
			y_cur[j] += tau * K_res[j];

		}
		res[i] = y_cur;
		prev_point = y_cur;
	}
	AddFileBatch(res, file, n); // последняя перезаписанная точка n-1
}

void write_file(string file, const vector<vector<double>>& GridFunc) {
	int num_vars = GridFunc[0].size();
	int n = GridFunc.size();

	ofstream out(file);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < num_vars; j++) {
			out << GridFunc[i][j] << " ";
		}
		out << endl;
	}
}


void write_file(string file, const vector<double>& vec) {
	int n = vec.size();

	ofstream out(file);
	for (int i = 0; i < n; i++) {
		out << vec[i] << endl;
	}
}



vector<vector<double>> Automated_Runge_Kutta(vector<double> y0,
	double T,
	double eps_tol,
	vector<double(*)(double, vector<double>)> syst, bool is_write_grid = 0,
	bool is_write_eps = 0, bool is_write_tau = 0) {
	vector<double> tau_n = { 1e-5 }; // задаём начальный шаг
	vector<vector<double>> res = { y0 };
	double cur_t = tau_n[0];
	int num_vars = syst.size(); // количество переменных
	const int p = 4; // порядок 
	int i = 0;
	//vector<double> k1_1(num_vars), k2_1(num_vars), k3_1(num_vars), k4_1(num_vars), K_1_res(num_vars);
	//vector<double> k1_2(num_vars), k2_2(num_vars), k3_2(num_vars), k4_2(num_vars), K_2_res(num_vars);
	//vector<double> y1_k1(num_vars), y1_k2(num_vars), y1_k3(num_vars),
	//    y2_k1(num_vars), y2_k2(num_vars), y2_k3(num_vars);
	vector<double> grid_tau, grid_half_tau;
	vector<double> grid, error, max_true_err;
	double true_x, true_y;
	grid.push_back(0);
	error.push_back(0);

	while (cur_t < T) {

		vector<double> y1, y2;
		grid_tau = { grid.back(), grid.back() + tau_n.back() }; // сетка с одним шагом
		grid_half_tau = {
			grid.back(),
			grid.back() + 0.5 * tau_n.back(),
			grid.back() + tau_n.back()
		}; // сетка из 2 шагов
		y1 = Runge_Kutta(res.back(), grid_tau, syst).back(); // решение, полученное на сетке из 1 шага
		y2 = Runge_Kutta(res.back(), grid_half_tau, syst).back();// решение, полученное на сетке из 2 шагов

		double estimate;
		estimate = inftyNorm(y1, y2) / (pow(2, p) - 1);
		if (estimate >= eps_tol) {
			tau_n.back() /= 2;

			cur_t = grid.back() + tau_n.back();
		}
		else if (estimate < 0.5 * eps_tol) {
			true_x = true_x_pendulum(cur_t);
			true_y = true_y_pendulum(cur_t);
			max_true_err.push_back(max(abs(true_x - y2[0]), abs(true_y - y2[1])));


			tau_n.push_back(2 * tau_n.back());
			//tau_n.push_back(tau_n.back());
			grid.push_back(cur_t);
			error.push_back(estimate);
			cur_t += tau_n.back();
			res.push_back(y2);

			i++;
			//PrintGridFunc(res);
		}
		else {
			true_x = true_x_pendulum(cur_t);
			true_y = true_y_pendulum(cur_t);
			max_true_err.push_back(max(abs(true_x - y2[0]), abs(true_y - y2[1])));


			grid.push_back(cur_t);
			error.push_back(estimate);
			tau_n.push_back(tau_n.back());
			//tau_n.push_back(0.5*tau_n.back());
			cur_t += tau_n.back();
			res.push_back(y2);
			i++;
			//PrintGridFunc(res);
		}


	}
	cout << *max_element(max_true_err.begin(), max_true_err.end());
	if (is_write_grid) {
		write_file("gridFile.txt", grid);
	}
	if (is_write_eps) {
		write_file("errorFile.txt", error);
	}
	if (is_write_tau) {
		write_file("tauFile.txt", tau_n);
	}
	return res;
}






vector<vector<double>> generate_2d_grid(double x0, double x1, double y0, double y1, int nx, int ny) {
	double hx = abs(x1 - x0) / nx, hy = abs(y1 - y0) / ny;
	double cx, cy;
	vector<vector<double>> res;
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


vector<vector<double>> generate_3d_grid(double x0, double x1, double y0, double y1, double z0, double z1, int nx, int ny, int nz) {
	double hx = abs(x1 - x0) / nx, hy = abs(y1 - y0) / ny, hz = abs(z1 - z0) / nz;
	double cx, cy, cz;
	vector<vector<double>> res;
	for (int i = 0; i <= ny; i++) {
		cx = x0;
		for (int j = 0; j <= nx; j++) {
			cz = z0;
			for (int k = 0; k <= nz; k++) {
				res.push_back({ cx, y0, cz });
				cz += hz;
			}

			cx += hx;

		}
		y0 += hy;

	}
	return res;
}


// center - центр квадрата, len - длина квадрата, k - кол-во точек в одну сторону
vector<vector<double>> GenerateStartPoints2D(const vector<double>& center, double len, double k) {
	int n = center.size(); // кол-во переменных
	vector<vector<double>> res;
	vector<double> x_i, newPoint(n);
	res = generate_2d_grid(center[0] - len / 2, center[0] + len / 2, center[1] - len / 2, center[1] + len / 2, k, k);
	return res;
}


// center - центр квадрата, len - длина квадрата, k - кол-во точек в одну сторону
vector<vector<double>> GenerateStartPoints3D(const vector<double>& center, double len, double k) {
	int n = center.size(); // кол-во переменных
	vector<vector<double>> res;
	vector<double> x_i, newPoint(n);
	res = generate_3d_grid(center[0] - len / 2, center[0] + len / 2, center[1] - len / 2, center[1] + len / 2,
		center[2] - len / 2, center[2] + len / 2, k, k, k);
	return res;
}


vector<vector<double>> DescribeArea(const vector<vector<double>>& start_points, function<vector<vector<double>>(vector<double> y0,
	const vector<double>& grid,
	vector<double(*)(double, vector<double>)> syst)> solver, vector<double(*)(double, vector<double>)> syst,
	vector<double> grid) {
	int n = start_points.size();
	vector<vector<double>> res;
	vector<vector<double>> temp;
	for (int i = 0; i < n; i++) {
		temp = solver(start_points[i], grid, syst);
		res.insert(res.end(), temp.begin(), temp.end());
	}
	return res;
}


double test4_1(double t,vector<double> vars)
{
	return vars[1];
}

//vector<double>& vars
double test4_2(double t, vector<double> vars)
{
	double M = 0.05;
	double delta = 0.01;
	double mu = 0.1;
	double nu = 0.1;
	double temp = (double)(M - 2 * delta * vars[1] - (mu + nu * sin(t)) * sin(vars[0]));
	return temp;
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
	//ProcessAitken(0.1, 10, { 2,0 }, { f1_system1_test , f2_system1_test }, 0.5, PredictionCorection);

	//WriteImplicitEuler(0.001, 10, { 0,0 }, { f1_system1_book , f2_system1_book });
	//vector<double> grid = GenerateUniformGrid(0, 1, 10);
	//vector<vector<double>> res = Euler_explicit({ 2,0 }, grid, { f1_system1_test, f2_system1_test });
	//PrintGridFunc(res);

	//res = SymmetricalScheme({ 2,0 }, grid, { f1_system1_test, f2_system1_test });
	//PrintGridFunc(res);

	//res = Runge_Kutta({ 2,0 }, grid, { f1_system1_test, f2_system1_test });
	//write_file("Runge_Kutta.txt", res);
	//vector<function<int(int)>> vec;
	//PrintGridFunc(res);


	// ===========================================================================
	// Маятник
	//vector<vector<double>> res;
	//vector<double> grid = GenerateUniformGrid(0, 1000, 1000);
	//res = ImplicitEuler(0.1, 10, { 1,0 }, { f1_pendulum, f2_pendulum });
	//write_file("ImplicitEulerPendulum.txt", res);
	//res = Euler_explicit({1,0},grid, { f1_pendulum, f2_pendulum });
	//write_file("ExplicitEulerPendulum.txt", res);
	//res = SymmetricalScheme({ 1,0 }, grid, { f1_pendulum, f2_pendulum });
	//write_file("SymmetricalSchemePendulum.txt", res);
	//res = Runge_Kutta({ 1,0 }, grid, { { f1_pendulum, f2_pendulum } });
	//write_file("RungeKuttaPendulum.txt", res);
	vector<vector<double>> res = PredictionCorection(0.01, 10, { 1,0 }, { f1_pendulum, f2_pendulum });
	write_file("PCpendulum.txt", res);

	// ===========================================================================
	// Фазовые портреты тестов из методички
	// Тест 1, точка (0,1)
	//vector<vector<double>> res,start_points;
	//start_points = GenerateStartPoints2D({ 0,1 }, 0.5, 20);
	//vector<double> grid = GenerateUniformGrid(0, 1, 10);
	//res = DescribeArea(start_points, Runge_Kutta, { f1_system1_book, f2_system1_book }, grid);
	//write_file("ExplicitEulerTest1P1.txt", res);
	//// Тест 1, точка (0,-1)
	//start_points = GenerateStartPoints2D({ 0,-1 }, 0.5, 10);
	//res = DescribeArea(start_points, Euler_explicit, { f1_system1_book, f2_system1_book }, grid);
	//write_file("ExplicitEulerTest1P2.txt", res);
	//// Тест 2 , точка (0,1)
	//start_points = GenerateStartPoints2D({ 0,1 }, 0.5, 10);
	//res = DescribeArea(start_points, Euler_explicit, { f1_system2_book, f2_system2_book }, grid);
	//write_file("ExplicitEulerTest2P1.txt", res);
	//// Тест 2 , точка (0,-1)
	//start_points = GenerateStartPoints2D({ 0,-1 }, 0.5, 20);
	//res = DescribeArea(start_points, Euler_explicit, { f1_system2_book, f2_system2_book }, grid);
	//write_file("ExplicitEulerTest2P2.txt", res);
	//// Тест 3
	//grid = GenerateUniformGrid(0, 1, 10);
	//start_points = GenerateStartPoints3D({ 1,1,1 }, 0.5, 10);
	//res = DescribeArea(start_points, SymmetricalScheme, { f1_system3_book, f2_system3_book, f3_system3_book }, grid);
	//write_file("SymmetricalSchemeTest3.txt", res);

	//============================================================================
	// Тест для функций с пакетом на маятнике
	//vector<vector<double>> res;
	//vector<double> grid = GenerateUniformGrid(0, 1000, 1e5);
	//BatchedEulerExplicit({ 1,0 }, grid, { f1_pendulum, f2_pendulum }, "pendulumExplicit.txt");
	//res = ImplicitEuler(0.01, 1000, { 1,0 }, { f1_pendulum, f2_pendulum });
	//write_file("pendulumImplicit.txt", res);
	//BatchedSymmetricalScheme({ 1,0 }, grid, { f1_pendulum, f2_pendulum }, "pendulumSymm.txt");
	//BatchedRungeKutta({ 1,0 }, grid, { f1_pendulum, f2_pendulum }, "pendulumRK.txt");
	//============================================================================
	// порядки сходимости
	/*for (double tau : {0.1,  0.01,  0.001, 0.0001}) {
		cout <<  tau << " & ";
		ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, Euler_explicit);
		cout << " & ";
		ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, ImplicitEuler);
		cout << " & ";
		ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, SymmetricalScheme);
		cout << " & ";
		ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, Runge_Kutta);
		cout << " & ";
	    ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, AdamsBashford);
	    cout << " & ";
	    ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, PredictionCorection);
	    cout << " \\\\ \\hline ";
	    cout << endl;

	}*/

	//cout << "Eitken" << endl;
	//for (double tau : {0.1, 0.01, 0.001, 0.0001}) {
	//	cout << tau << endl;
	//	ProcessAitken(tau, 10, { 2,0 }, { f1_variant , f2_variant }, 0.5, ImplicitEuler);
	//	cout << " ";
	//	ProcessAitken(tau, 10, { 2,0 }, { f1_variant , f2_variant }, 0.5, AdamsBashford);
	//	cout << " ";
	//	ProcessAitken(tau, 10, { 2,0 }, { f1_variant , f2_variant }, 0.5, PredictionCorection);
	//	cout << endl;


	//}
	//cout << "AnalitResh" << endl;
	//for (double tau : {0.1, 0.01, 0.001, 0.0001}) {
	//	cout << tau << endl;
	//	PorydokAnal(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, ImplicitEuler, sol_pendulum);
	//	cout << " ";
	//	PorydokAnal(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, AdamsBashford, sol_pendulum);
	//	cout << " ";
	//	PorydokAnal(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, PredictionCorection, sol_pendulum);
	//	cout << endl;


	//}
	//for (double tau : {0.1, 0.01, 0.001, 0.0001}) {
	//	cout << tau << " & ";
	//	//ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, Euler_explicit);
	//	//cout << " & ";
	//	GetOrder(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, ImplicitEuler,{true_x_pendulum, true_y_pendulum});
	//	cout << " & ";
	//	/*ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, SymmetricalScheme);
	//	cout << " & ";
	//	ProcessAitken(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, Runge_Kutta);
	//	cout << " & ";*/
	//	GetOrder(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, AdamsBashford,  { true_x_pendulum, true_y_pendulum });
	//	cout << " & ";
	//	GetOrder(tau, 10, { 2,0 }, { f1_pendulum , f2_pendulum }, 0.5, PredictionCorection,  { true_x_pendulum, true_y_pendulum });
	//	cout << " \\\\ \\hline ";
	//	cout << endl;

	//}
	//cout << endl;
	//ProcessAitken(0.01, 10, { 1,0 }, { f1_pendulum, f2_pendulum }, 0.5, Euler_explicit);
	//cout << endl;
	//ProcessAitken(0.01, 10, { 1,0 }, { f1_pendulum, f2_pendulum }, 0.5, ImplicitEuler);
	//cout << endl;
	//ProcessAitken(0.001, 10, { 1,0 }, { f1_pendulum, f2_pendulum }, 0.5, SymmetricalScheme);
	//cout << endl;
	//ProcessAitken(0.01, 10, { 1,0 }, { f1_pendulum, f2_pendulum }, 0.5, Runge_Kutta);
	// =============================================================================================
	// Автоматический шаг. Маятник

	//vector<vector<double>> res = Automated_Runge_Kutta({ 0.1,0.1 }, 100, 1e-6, { test4_1, test4_2 }, 1, 1, 1);
	//write_file("AutomatedStep2.txt", res);
	//vector<double> grid = GenerateUniformGrid(0, 10, 1e6);
	//write_file("gridFile2.txt", grid);
	//vector<vector<double>> res = Runge_Kutta({ 1,0 }, grid, { f1_pendulum, f2_pendulum });
	//write_file("test123.txt", res);
	// =============================================================================================
	// Система по вариантам
	//vector<double> grid = GenerateUniformGrid(0, 200, 1000);
	//vector<vector<double>> res = Runge_Kutta({ 0.1, 0.1 }, grid, { f1_variant, f2_variant });
	//write_file("variant_syst.txt", res);
	//vector<vector<double>> start_points = GenerateStartPoints2D({ 0,0 }, 0.5, 20);
	//grid = GenerateUniformGrid(0, 1, 20);;
	//res = DescribeArea(start_points, Runge_Kutta, { f1_variant, f2_variant }, grid);
	//write_file("variant_syst_SP.txt", res);
}
