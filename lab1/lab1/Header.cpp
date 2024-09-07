#include "Header.h"


void triangularize(vector<vector<double>>& matrix) {
	int n = matrix.size();

	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			double temp = matrix[i][j];
			for (int k = 0; k < n+1; ++k) {
				matrix[i][k] -=  matrix[j][k] * temp / matrix[j][j];
			}
		}
	}
}


void print_matrix(const vector<vector<double>>& matrix) {
	int n = matrix.size();
	cout << n << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n+1; ++j) {
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
			for (int j = 0; j < n+1; ++j) {
				out << matrix[i][j] << " ";
			}
			out << endl;
		}
		out.close();
	}
}



vector<vector<double>> readf_matrix(string file_matrix, string file_vector) {
	vector<double> b = readf_vector(file_vector);
	ifstream in;
	in.open(file_matrix);
	int n;
	vector<vector<double>> matrix;

	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			matrix.emplace_back(vector<double>(n+1));
			for (int j = 0; j < n; ++j) {
				in >> matrix[i][j];
			}
			matrix[i][n] = b[i];
		}
		in.close();
	}
	return matrix;
}


int get_main_element_row(const vector<vector<double>>& matrix, int var_row) {
	double max = matrix[var_row][var_row];
	double max_row = var_row;
	int n = matrix.size();
	for (int i = var_row + 1; i < n; ++i) {
		if (matrix[i][var_row] > max) {
			max = matrix[i][var_row];
			max_row = i;
		}
	}
	return max_row;
}


vector<double> readf_vector(string file) {
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


vector<double> Gauss_method(vector<vector<double>>& matrix){
	permutate_rows(matrix);
	print_matrix(matrix);
	triangularize(matrix);
	print_matrix(matrix);
	int n = matrix.size();
	vector<double> x(n);
	double sum;
	for (int i = n - 1; i > -1; --i) {
		sum = 0;
		for (int j = i + 1; j < n; ++j) {
			sum += matrix[i][j] * x[j];
		}
		x[i] = (matrix[i][n] - sum) / matrix[i][i];
		//cout << x[i];
	}
	return x;
}


vector<double> check_ans(const vector<vector<double>>& matrix, const vector<double>& x) {
	int n = matrix.size();
	vector<double> diff(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			diff[i] += matrix[i][j] * x[j];
		}
		diff[i] -= matrix[i][n];
	}
	return diff;
}