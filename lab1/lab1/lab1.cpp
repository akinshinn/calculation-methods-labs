#include "Gauss.h"

void main()
{
	vector<vector<double>> a = readf_matrix_double("data\\matrix5.txt");
	print_matrix(a);
	vector<double> x = Gauss_method(a);
	for (int i = 0; i < x.size(); ++i) cout << x[i] << " ";
	write_file("out1.txt", x);
	DoubleTest("LinearSystem.txt");
}


