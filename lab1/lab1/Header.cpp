#include "Header.h"


void get_U(vector<vector<double>>* matrix) {
	int n = (*matrix).size();

	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			for (int k = 0; k < n; ++k) {
				(*matrix)[i][k] = (*matrix)[i][k] - (*matrix)[j][k] * (*matrix)[i][j] / (*matrix)[j][j];
			}
		}
	}
}


void print_matrix(const vector<vector<double>>& matrix) {
	int n = matrix.size();
	cout << n << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


vector<vector<double>> read_matrix() {
	int n; 
	cin >> n;
	vector<vector<double>> matrix;
	
	for (int i = 0; i < n; ++i) {
		matrix.emplace_back(vector<double>(n));
		for (int j = 0; j < n; ++j) {
			cin >> matrix[i][j];
		}
	}
	return matrix;
}


void write_file(const vector<vector<double>>& matrix, string file) {
	ofstream out;
	out.open(file);
	if (out.is_open()) {
		int n = matrix.size();
		out << n << endl;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				out << matrix[i][j] << " ";
			}
			out << endl;
		}
		out.close();
	}
}



vector<vector<double>> readf_matrix(string file) {
	ifstream in;
	in.open(file);
	int n;
	vector<vector<double>> matrix;

	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			matrix.emplace_back(vector<double>(n));
			for (int j = 0; j < n; ++j) {
				in >> matrix[i][j];
			}
		}
		in.close();
	}
	return matrix;
}