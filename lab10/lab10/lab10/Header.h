#pragma once
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
class Matrix {
public:
	vector<vector<double>> elements;
	int size;
	Matrix(int n) {
		size = n;
		for (int i = 0; i < n; i++) {
			elements.push_back(vector<double>(n, 0));
		}
	}
	Matrix(const vector<vector<double>>& matrix) {
		size = matrix.size();
		elements = matrix;
	}
	vector<double>& operator[](int i) {
		return elements[i];
	}
	vector<double> operator[](int i) const {
		return elements[i];
	}

	Matrix& operator*(double x) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				elements[i][j] *= x;
			}
		}
		return *this;
	}
	Matrix& operator=(Matrix m) {
		size = m.size;
		elements = m.elements;
		return *this;
	}
	Matrix operator*(const Matrix& m) const {
		return product(*this, m);
	}
	Matrix operator+(const Matrix& m) const {
		return sum(*this, m);
	}
	Matrix operator-(const Matrix& m) const {
		return diff(*this, m);
	}
	void print() const {
		cout << endl;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				cout << elements[i][j] << " " << setprecision(20);
			}
			cout << endl;
		}
	}
	friend Matrix product(const Matrix& m1, const Matrix& m2);
	friend Matrix sum(const Matrix& m1, const Matrix& m2);
	friend Matrix diff(const Matrix& m1, const Matrix& m2);
};


Matrix product(const Matrix& m1, const Matrix& m2) {
	int n = m1.size;
	Matrix res(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0;
			for (int k = 0; k < n; k++) {
				sum += m1[i][k] * m2[k][j];
			}
			res[i][j] = sum;
		}
	}
	return res;
}


Matrix sum(const Matrix& m1, const Matrix& m2) {
	int n = m1.size;
	Matrix res(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res[i][j] = m1[i][j] + m2[i][j];
		}
	}
	return res;
}

Matrix diff(const Matrix& m1, const Matrix& m2) {
	int n = m1.size;
	Matrix res(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res[i][j] = m1[i][j] - m2[i][j];
		}
	}
	return res;
}





Matrix eye(int n) {
	Matrix res(n);
	for (int i = 0; i < n; i++) {
		res[i][i] = 1;
	}
	return res;
}


void triangularize(Matrix& matrix, vector<double>& b) {
	int n = matrix.size;

	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			double temp = matrix[i][j];
			for (int k = 0; k < n; ++k) {
				matrix[i][k] -= matrix[j][k] * temp / matrix[j][j];
			}
			b[i] -= b[j] * temp / matrix[j][j];
		}
	}
}


bool check_matrix(const Matrix& triang_matrix) {
	double num = 1;
	for (int i = 0; i < triang_matrix.size; ++i) {
		num *= triang_matrix[i][i];
	}
	return (abs(num) > 1e-15);
}

vector<double> Gauss_method(Matrix matrix, vector<double> b) {
	triangularize(matrix, b);
	//matrix.print();
	int n = matrix.size;
	vector<double> x(n);
	if (!check_matrix(matrix)) {
		cout << "det A = 0" << endl;
		//matrix.print();
		return x;
	}

	double sum;
	for (int i = n - 1; i > -1; --i) {
		sum = 0;
		for (int j = i + 1; j < n; ++j) {
			sum += matrix[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / matrix[i][i];
	}
	return x;
}