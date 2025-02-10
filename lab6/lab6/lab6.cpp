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
using namespace std;

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


vector<double>NewtonMethod(vector<double> y0, vector<double(&f)(vector<double> var)>) {


    return f;
}

int main()
{
    
}
