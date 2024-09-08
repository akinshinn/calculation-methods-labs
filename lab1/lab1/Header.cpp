#include "Header.h"


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
vector<float> readf_vector_float(string file) {
	ifstream in;
	in.open(file);
	int n;
	vector<float> b;
	float temp;
	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			in >> temp;
			b.emplace_back(temp);
		}
		in.close();
	}
	return b;
}

vector<vector<float>> readf_matrix_float(string file_matrix, string file_vector) {
	vector<float> b = readf_vector_float(file_vector);
	ifstream in;
	in.open(file_matrix);
	int n;
	vector<vector<float>> matrix;

	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			matrix.emplace_back(vector<float>(n + 1));
			for (int j = 0; j < n; ++j) {
				in >> matrix[i][j];
			}
			matrix[i][n] = b[i];
		}
		in.close();
	}
	return matrix;
}









vector<double> readf_vector_double(string file) {
	ifstream in;
	in.open(file);
	int n;
	vector<double> b;
	double temp;
	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			in >> temp;
			b.emplace_back(temp);
		}
		in.close();
	}
	return b;
}




vector<vector<double>> readf_matrix_double(string file_matrix, string file_vector) {
	vector<double> b = readf_vector_double(file_vector);
	ifstream in;
	in.open(file_matrix);
	int n;
	vector<vector<double>> matrix;

	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			matrix.emplace_back(vector<double>(n + 1));
			for (int j = 0; j < n; ++j) {
				in >> matrix[i][j];
			}
			matrix[i][n] = b[i];
		}
		in.close();
	}
	return matrix;
}


