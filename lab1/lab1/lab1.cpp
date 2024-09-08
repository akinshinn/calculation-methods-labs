#include "Header.h"

int main()
{
	vector<vector<float>> a = readf_matrix_float("data\\matrix.txt", "data\\vector.txt");
	//readf_vector_float("data\\vector.txt");
	print_matrix(a);
	vector<vector<float>> b = a;
	vector<float> x = Gauss_method(a);
	for (int i = 0; i < x.size(); ++i) cout << x[i] << " ";
	cout << endl;
	vector<float> diff = check_ans(a, x);
	for (int i = 0; i < diff.size(); ++i) cout << diff[i] << " ";
	vector<vector<float>> inv = Gauss_inverse(b);
	print_matrix(inv);
	vector<vector<float>> prod = matrix_prod(b, inv);
	print_matrix(prod);
	write_file(b, x, inv, "data\\out.txt");
}
