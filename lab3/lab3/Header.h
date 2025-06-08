#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//#include "QR.cpp"
//#pragma once
//#include <vector>
//#include <iostream>
//#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>

using namespace std;


pair<vector<double>, int> QRDoubleDecomposition(int size, vector<vector<double>>& Matrix, vector<vector<double>>& TMatrix) {
    // ������ ���
    double c, s;
    double tmp1, tmp2, tmp3, tmp4;
    for (int i = 0; i < size - 1; i++) {
        // i - index ���������� ������� ����� ���������, ���� i = 2 ��������� x2 �� ���� ���������
        for (int j = i + 1; j < size; j++) {
            // j - index ��������� � ������� ����� �������� x_i
            if (Matrix[j][i] == 0) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);

            for (int k = 0; k < size + 1; k++) {
                // size + 1 ������ ��� Matrix ������� �� n+1 ������� Matrix = A | b
                tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
                tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
                Matrix[i][k] = tmp1;
                Matrix[j][k] = tmp2;
            }


            for (int k = 0; k < size; k++) {
                tmp3 = c * TMatrix[i][k] + s * TMatrix[j][k];
                tmp4 = -s * TMatrix[i][k] + c * TMatrix[j][k];
                TMatrix[i][k] = tmp3;
                TMatrix[j][k] = tmp4;
            }
        }
    }
    // ��������
    double product = 1;
    for (int i = 0; i < size; i++) {
        product *= Matrix[i][i];
    }

    vector<double> ans(size);
    int flag = 1;
    if (product != 0) {
        // �������� ���
        int n = size - 1;
        double suma = 0;
        for (int k = 0; k < size; k++) {
            suma = 0;
            for (int i = n - k + 1; i < size; i++) {
                suma += ans[i] * Matrix[n - k][i];
            }
            ans[n - k] = (Matrix[n - k][n + 1] - suma) / Matrix[n - k][n - k];
        }
    }
    else {
        flag = 0;
    }

    return { ans , flag };
}
