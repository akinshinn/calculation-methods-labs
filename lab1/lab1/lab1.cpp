#include "Header.h"

int main()
{
	vector<vector<double>> a = readf_matrix("data\\matrix.txt", "data\\vector.txt");
	print_matrix(a);
	vector<double> x = Gauss_method(a);
	for (int i = 0; i < x.size(); ++i) cout << x[i] << " ";
	cout << endl;
	vector<double> diff = check_ans(a, x);
	for (int i = 0; i < diff.size(); ++i) cout << diff[i] << " ";
	//write_file(a, "data\\out.txt");
}
