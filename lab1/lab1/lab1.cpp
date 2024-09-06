#include "Header.h"

int main()
{
	vector<vector<double>> a = readf_matrix("data\\in.txt");
	print_matrix(a);
	get_U(&a);
	print_matrix(a);
	write_file(a, "data\\out.txt");
	
}
