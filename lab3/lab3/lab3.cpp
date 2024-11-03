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
const double epsilon = 0.01;
void DisplayMatrix(vector<vector<double>> Matrix) {
    for (auto& row : Matrix) {
        for (auto el : row) {
            cout << el << " ";
        }
        cout << endl;
    }
}


void DisplayVector(vector<double> vec) {
    for (auto el : vec) {
        cout << el << " ";
    }
    cout << endl;
}


vector<vector<double>> ProductMatrix(vector<vector<double>> mat1, vector<vector<double>> mat2) {
    vector<vector<double>> res = {};
    vector<double> str{};
    int n = mat1.size();
    double sum = 0;
    for (int i = 0; i < n; i++) {
        str = {};
        for (int j = 0; j < n; j++) {
            sum = 0;
            for (int k = 0; k < n; k++) {
                sum += mat1[i][k] * mat2[k][j];
            }
            str.push_back(sum);
        }
        res.push_back(str);
    }
    return res;
}


bool StopCriteria(vector<vector<double>> Matrix) {
    double max_el = -numeric_limits<double>::max();
    for (int i = 0; i < Matrix.size(); i++) {
        for (int j = 0; j < Matrix.size(); j++) {
            if (i > j) {
                max_el = max(max_el, abs(Matrix[i][j]));
            }
        }
    }
    return max_el < epsilon;
}

pair<vector<vector<double>>, vector<vector<double>>> QRDecomposition(vector<vector<double>> Matrix, int& prod) {

    double c, s;
    int size = Matrix.size();
    double tmp1, tmp2, tmp3, tmp4;
    vector<vector<double>> TMatrix, QMatrix;
    vector<double> test;
    for (int i = 0; i < Matrix.size(); i++){
        test = {};
        for (int j = 0; j < Matrix.size(); j++) {
            if (i == j) { test.push_back(1); }
            else { test.push_back(0); }

        }
        TMatrix.push_back(test);
    }

    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        for (int j = i + 1; j < size; j++) {
            // j - index уравнения в котором будем устанять x_i
            if (Matrix[j][i] == 0) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            prod += 4;

            for (int k = i; k < size; k++) {
                tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
                tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
                prod += 4;
                Matrix[i][k] = tmp1;
                Matrix[j][k] = tmp2;
            }

            for (int k = 0; k < size; k++) {
                tmp3 = c * TMatrix[i][k] + s * TMatrix[j][k];
                tmp4 = -s * TMatrix[i][k] + c * TMatrix[j][k];
                prod += 4;
                TMatrix[i][k] = tmp3;
                TMatrix[j][k] = tmp4;
            }
        }
    }
    
    for (int i = 0; i < Matrix.size(); i++) {
        test = {};
        for (int j = 0; j < Matrix.size(); j++) {
            test.push_back(TMatrix[j][i]);

        }
        QMatrix.push_back(test);
    }

    return { QMatrix, Matrix };
}


pair<vector<vector<double>>, vector<vector<double>>> QRDecompositionForShift(vector<vector<double>> Matrix,int size, int& prod) {

    double c, s;
    double tmp1, tmp2, tmp3, tmp4;
    vector<vector<double>> TMatrix, QMatrix;
    vector<double> test;
    for (int i = 0; i < size; i++) {
        test = {};
        for (int j = 0; j < size; j++) {
            if (i == j) { test.push_back(1); }
            else { test.push_back(0); }

        }
        TMatrix.push_back(test);
    }

    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        for (int j = i + 1; j < size; j++) {
            // j - index уравнения в котором будем устанять x_i
            if (Matrix[j][i] == 0) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            prod += 4;

            for (int k = i; k < size; k++) {
                tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
                tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
                prod += 4;
                Matrix[i][k] = tmp1;
                Matrix[j][k] = tmp2;
            }

            for (int k = 0; k < size; k++) {
                tmp3 = c * TMatrix[i][k] + s * TMatrix[j][k];
                tmp4 = -s * TMatrix[i][k] + c * TMatrix[j][k];
                prod += 4;
                TMatrix[i][k] = tmp3;
                TMatrix[j][k] = tmp4;
            }
        }
    }

    for (int i = 0; i < size; i++) {
        test = {};
        for (int j = 0; j < size; j++) {
            test.push_back(TMatrix[j][i]);

        }
        QMatrix.push_back(test);
    }

    return { QMatrix, Matrix };
}


pair<vector<vector<double>>, vector<vector<double>>> QRDecompositionForHesenberg(vector<vector<double>> Matrix, int& prod) {
    double c, s;
    int size = Matrix.size();
    double tmp1, tmp2, tmp3, tmp4;
    vector<vector<double>> TMatrix, QMatrix;
    vector<double> test;
    for (int i = 0; i < Matrix.size(); i++) {
        test = {};
        for (int j = 0; j < Matrix.size(); j++) {
            if (i == j) { test.push_back(1); }
            else { test.push_back(0); }

        }
        TMatrix.push_back(test);
    }

    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        int j = i + 1;
        // j - index уравнения в котором будем устанять x_i
        c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
        s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
        prod += 4;

        for (int k = i; k < size; k++) {
            tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
            tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
            prod += 4;
            Matrix[i][k] = tmp1;
            Matrix[j][k] = tmp2;
        }

        for (int k = 0; k < size; k++) {
            tmp3 = c * TMatrix[i][k] + s * TMatrix[j][k];
            tmp4 = -s * TMatrix[i][k] + c * TMatrix[j][k];
            prod += 4;
            TMatrix[i][k] = tmp3;
            TMatrix[j][k] = tmp4;
        }
        
    }

    for (int i = 0; i < Matrix.size(); i++) {
        test = {};
        for (int j = 0; j < Matrix.size(); j++) {
            test.push_back(TMatrix[j][i]);

        }
        QMatrix.push_back(test);
    }

    return { QMatrix, Matrix };
}


pair<vector<vector<double>>, vector<vector<double>>> QRDecompositionForHesenbergAndShift(vector<vector<double>> Matrix,int size, int& prod) {
    double c, s;
    double tmp1, tmp2, tmp3, tmp4;
    vector<vector<double>> TMatrix, QMatrix;
    vector<double> test;
    for (int i = 0; i < size; i++) {
        test = {};
        for (int j = 0; j < size; j++) {
            if (i == j) { test.push_back(1); }
            else { test.push_back(0); }

        }
        TMatrix.push_back(test);
    }

    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        int j = i + 1;
        // j - index уравнения в котором будем устанять x_i
        c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
        s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
        prod += 4;

        for (int k = i; k < size; k++) {
            tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
            tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
            prod += 4;
            Matrix[i][k] = tmp1;
            Matrix[j][k] = tmp2;
        }

        for (int k = 0; k < size; k++) {
            tmp3 = c * TMatrix[i][k] + s * TMatrix[j][k];
            tmp4 = -s * TMatrix[i][k] + c * TMatrix[j][k];
            prod += 4;
            TMatrix[i][k] = tmp3;
            TMatrix[j][k] = tmp4;
        }

    }

    for (int i = 0; i < size; i++) {
        test = {};
        for (int j = 0; j < size; j++) {
            test.push_back(TMatrix[j][i]);

        }
        QMatrix.push_back(test);
    }

    return { QMatrix, Matrix };
}


void DataRead(vector<vector<double>>& Matrix, string filename) {
    ifstream file(filename);
    string line;
    int intSize = 0;
    getline(file, line);
    stringstream str(line);
    str >> intSize;
    Matrix.resize(intSize, vector<double>(intSize, 0));
    for (auto& row : Matrix) {
        row.resize(intSize, 0);
    }

    for (int Curstr = 0; Curstr < intSize; Curstr++) {
        getline(file, line);
        stringstream str(line);
        for (int j = 0; j < intSize; j++) {
            str >> Matrix[Curstr][j];
        }
    }

    file.close();
}


vector<double> EigenValuesQR(vector<vector<double>> Matrix, int& prod) {
    pair<vector<vector<double>>, vector<vector<double>>> QR;
    int iter = 0,n;
    n = Matrix.size();
    while (!StopCriteria(Matrix)) {
        QR = QRDecomposition(Matrix, prod);
        for (int i = 0; i < Matrix.size(); i++) {
            for (int j = 0; j < Matrix.size(); j++) {
                Matrix[i][j] = 0;
                for (int k = i; k < Matrix.size(); k++) {
                    Matrix[i][j] += QR.second[i][k] * QR.first[k][j];
                    prod++;
                }
            }
        }
        iter++;
        //DisplayMatrix(Matrix);
        //cout << endl;
    }


    cout << "Iterations: " << iter << endl;
    // можно считать по формуле, но это плохо тк если я код поменяю надо формулы пересчитывать
    //cout << "Products: " << iter * n * (9 * n * n - 3 * n - 4) / 2 << endl;
    cout << "Products: " << prod << endl;
    vector<double> EigenValues{};
    for (int i = 0; i < Matrix.size(); i++) {
        EigenValues.push_back(Matrix[i][i]);

    }
    return EigenValues;
}


vector<vector<double>> HesenbergDecomposition(vector<vector<double>> Matrix, int& prod) {
    
    double alpha, beta, koren,tmp1,tmp2;
    for (int k = 1; k < Matrix.size() - 1; k++) {
        for (int l = k+1; l < Matrix.size(); l++) {
            koren = sqrt(Matrix[l][k - 1] * Matrix[l][k - 1] + Matrix[k][k - 1] * Matrix[k][k - 1]);
            alpha = Matrix[k][k - 1] / koren;
            beta = Matrix[l][k - 1] / koren;
            prod += 4;
            for (int index = k-1; index < Matrix.size(); index++) {
                tmp1 = alpha * Matrix[k][index] + beta * Matrix[l][index];
                tmp2 = alpha * Matrix[l][index] - beta * Matrix[k][index];
                Matrix[k][index] = tmp1;
                Matrix[l][index] = tmp2;
                prod += 4;
            }
            for (int index = 0; index < Matrix.size(); index++) {
                tmp1 = alpha * Matrix[index][k] + beta * Matrix[index][l];
                tmp2 = alpha * Matrix[index][l] - beta * Matrix[index][k];
                Matrix[index][k] = tmp1;
                Matrix[index][l] = tmp2;
                prod += 4;
            }
        }
    }
    return Matrix;

}


vector<double> EigenValuesHesenberg(vector<vector<double>> Matrix,int& prod) {
    vector<double> EigenValues = {};
    Matrix = HesenbergDecomposition(Matrix, prod);
    pair<vector<vector<double>>, vector<vector<double>>> QR;
    int iter = 0, n;
    n = Matrix.size();
    while (!StopCriteria(Matrix)) {
        QR = QRDecompositionForHesenberg(Matrix, prod);
        for (int i = 0; i < Matrix.size(); i++) {
            for (int j = 0; j < Matrix.size(); j++) {
                Matrix[i][j] = 0;
                for (int k = i; k < Matrix.size(); k++) {
                    Matrix[i][j] += QR.second[i][k] * QR.first[k][j];
                    prod++;
                }
            }
        }
        iter++;
        //DisplayMatrix(Matrix);
        //cout << endl;
    }


    cout << "Iterations: " << iter << endl;
    // можно считать по формуле, но это плохо тк если я код поменяю надо формулы пересчитывать
    //cout << "Products: " << iter * n * (9 * n * n - 3 * n - 4) / 2 << endl;
    cout << "Products: " << prod << endl;
    for (int i = 0; i < Matrix.size(); i++) {
        EigenValues.push_back(Matrix[i][i]);
    }
    return EigenValues;
}

bool StopCriteriaShift(vector<vector<double>> Matrix, int curr) {
    double max_el = -numeric_limits<double>::max();
    int size = Matrix.size();
    int index = size - curr - 1;
    for (int j = 0; j < index; j++) {
        max_el = max(max_el, abs(Matrix[index][j]));
    }
    return max_el < epsilon;
}

vector<double> EigenValuesShift(vector<vector<double>> Matrix, int& prod) {
    pair<vector<vector<double>>, vector<vector<double>>> QR;
    vector<double> res{};
    int iter = 0, size = Matrix.size(), curr =0 ;
    // curr - количество найденных собственных значений 
    while (curr < size) {
        // сдвигаем A - sigma E, e; уже найденные собственные значения не трогаем
        double sigma = Matrix[size - curr - 1][size - curr - 1];
        for (int i = 0; i < size - curr; i++) {
            Matrix[i][i] -= sigma;
        }

        // ищем с помощью обычного QR разложения собственные значения.
        while (!StopCriteriaShift(Matrix, curr)) {
            QR = QRDecompositionForShift(Matrix, size - curr, prod);
            for (int i = 0; i < size-curr; i++) {
                for (int j = 0; j < size - curr; j++) {
                    Matrix[i][j] = 0;
                    for (int k = i; k < size - curr; k++) {
                        Matrix[i][j] += QR.second[i][k] * QR.first[k][j];
                        prod++;
                    }
                }
            }
            /*DisplayMatrix(Matrix);
            cout << endl;*/
            iter++;
        }
        /*DisplayMatrix(Matrix);
        cout << endl;*/
        //сдвигаем обратно
        for (int i = 0; i < size - curr; i++) {
            Matrix[i][i] += sigma;
        }
        curr++;
    }
    // собираем все собственные значения 
    for (int i = 0; i < size; i++) {
        res.push_back(Matrix[i][i]);
    }

    cout << "Iterations: " << iter << endl;
    cout << "Products: " << prod << endl;

    return res;
}


vector<double> EigenValuesHesenbergAndShift(vector<vector<double>> Matrix, int& prod) {
    vector<double> EigenValues = {};
    Matrix = HesenbergDecomposition(Matrix, prod);
    pair<vector<vector<double>>, vector<vector<double>>> QR;
    vector<double> res{};
    int iter = 0, size = Matrix.size(), curr = 0;
    // curr - количество найденных собственных значений 
    while (curr < size) {
        // сдвигаем A - sigma E, e; уже найденные собственные значения не трогаем
        double sigma = Matrix[size - curr - 1][size - curr - 1];
        for (int i = 0; i < size - curr; i++) {
            Matrix[i][i] -= sigma;
        }

        // ищем с помощью обычного QR разложения собственные значения.
        while (!StopCriteriaShift(Matrix, curr)) {
            QR = QRDecompositionForHesenbergAndShift(Matrix, size - curr, prod);
            for (int i = 0; i < size - curr; i++) {
                for (int j = 0; j < size - curr; j++) {
                    Matrix[i][j] = 0;
                    for (int k = i; k < size - curr; k++) {
                        Matrix[i][j] += QR.second[i][k] * QR.first[k][j];
                        prod++;
                    }
                }
            }
            /*DisplayMatrix(Matrix);
            cout << endl;*/
            iter++;
        }
        /*DisplayMatrix(Matrix);
        cout << endl;*/
        //сдвигаем обратно
        for (int i = 0; i < size - curr; i++) {
            Matrix[i][i] += sigma;
        }
        curr++;
    }
    // собираем все собственные значения 
    for (int i = 0; i < size; i++) {
        res.push_back(Matrix[i][i]);
    }

    cout << "Iterations: " << iter << endl;
    cout << "Products: " << prod << endl;

    return res;
}
int main() {

    vector<vector<double>> Matrix;
    DataRead(Matrix, "Matrix.txt");
    int prodVanilaQr = 0, prodQrWithHesenberg=0, prodQrWithShift=0, prodQrWithHesenbergAndShift=0;

    cout << "Vanila QR decompositon" << endl;
    DisplayVector(EigenValuesQR(Matrix, prodVanilaQr));

    cout << endl;
    cout << "QR with Hesenberg" << endl;
    DisplayVector(EigenValuesHesenberg(Matrix, prodQrWithHesenberg));

    cout << endl;
    cout << "QR with Shift" << endl;
    DisplayVector(EigenValuesShift(Matrix, prodQrWithShift));

    cout << endl;
    cout << "QR with Hesenberg And Shift" << endl;
    DisplayVector(EigenValuesHesenbergAndShift(Matrix, prodQrWithHesenbergAndShift));
}