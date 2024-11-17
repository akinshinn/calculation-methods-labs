#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>

using namespace std;
const double epsilon = 1e-4;
const int MAX_ITERATIONS = 1000;
const double DELTA = 1e-2;
double hessenberg = 0;

void DisplayMatrix(vector<vector<double>> Matrix) {
    for (auto& row : Matrix) {
        for (auto el : row) {
            cout << el << " ";
        }
        cout << endl;
    }
}


double NormInfVec(const vector<double>& vec) {
    double Maxel = -1;
    for (int i = 0; i < vec.size(); i++) {
        Maxel = max(Maxel, abs(vec[i]));
    }
    return Maxel;
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
    double max_el = -1;
    int n = Matrix.size();
    for (int i = 0; i < Matrix.size()-1; i++) {
        for (int j = 0; j < Matrix.size()-1; j++) {
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


pair<vector<vector<double>>, vector<vector<double>>> QRDecompositionForShift(vector<vector<double>> Matrix, int size, int& prod) {

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
    //cout << "Hessengberg products: " << hessenberg << endl;
    return { QMatrix, Matrix };
}


pair<vector<vector<double>>, vector<vector<double>>> QRDecompositionForHesenbergAndShift(vector<vector<double>> Matrix, int size, int& prod) {
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
    int iter = 0, n;
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
        cout << "iteration " << iter << " " << prod << endl;
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

    double alpha, beta, koren, tmp1, tmp2;
    for (int k = 1; k < Matrix.size() - 1; k++) {
        for (int l = k + 1; l < Matrix.size(); l++) {
            koren = sqrt(Matrix[l][k - 1] * Matrix[l][k - 1] + Matrix[k][k - 1] * Matrix[k][k - 1]);
            alpha = Matrix[k][k - 1] / koren;
            beta = Matrix[l][k - 1] / koren;
            prod += 4;
            hessenberg += 4;

            for (int index = k - 1; index < Matrix.size(); index++) {
                tmp1 = alpha * Matrix[k][index] + beta * Matrix[l][index];
                tmp2 = alpha * Matrix[l][index] - beta * Matrix[k][index];
                Matrix[k][index] = tmp1;
                Matrix[l][index] = tmp2;
                prod += 4;
                hessenberg += 4;

            }
            for (int index = 0; index < Matrix.size(); index++) {
                tmp1 = alpha * Matrix[index][k] + beta * Matrix[index][l];
                tmp2 = alpha * Matrix[index][l] - beta * Matrix[index][k];
                Matrix[index][k] = tmp1;
                Matrix[index][l] = tmp2;
                prod += 4;
                hessenberg += 4;

            }
        }
    }
    //cout << "products " << hessenberg << endl;
    return Matrix;

}


vector<double> EigenValuesHesenberg(vector<vector<double>> Matrix, int& prod) {
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
                    //prod++; 
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
    int iter = 0, size = Matrix.size(), curr = 0;
    // curr - количество найденных собственных значений 
    while (curr < size) {
        // сдвигаем A - sigma E, e; уже найденные собственные значения не трогаем
        double sigma = Matrix[size - curr - 1][size - curr - 1];
        for (int i = 0; i < size - curr; i++) {
            Matrix[i][i] -= sigma;
        }

        // ищем с помощью обычного QR разложения собственные значения.
        //while (!StopCriteriaShift(Matrix, curr)) {
        while (!StopCriteriaShift(Matrix, curr)) {
            QR = QRDecompositionForShift(Matrix, size - curr, prod);
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

vector<double> matrix_prod_vec(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    vector<double> res(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int k = 0; k < n; ++k) {
            sum += A[i][k] * b[k];
        }
        res[i] = sum;
    }
    return res;
}


double eigen_check(vector<vector<double>> ModifiedA, vector<double> vect) {
    vector<double> vec1 = matrix_prod_vec(ModifiedA, vect);
    //DisplayMatrix(ModifiedA);
    //cout << NormInfVec(vec_diff(vec1, vec2));
    return NormInfVec(vec1);
}

double NormSquareVec(const vector<double>& vec) {
    double suma = 0;
    for (int i = 0; i < vec.size(); i++) {
        suma += vec[i] * vec[i];
    }
    return sqrt(suma);
}



vector<vector<double>> transpose(const vector<vector<double>>& matrix) {
    vector<vector<double>> res = matrix;
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            res[i][j] = matrix[j][i];
            res[j][i] = matrix[i][j];
        }
    }
    return res;
}

vector<double> QR_final_step(const vector<vector<double>>& Q_transpose, const vector<vector<double>>& R, vector<double> b) {
    vector<double> res(R.size());

    b = matrix_prod_vec(Q_transpose, b);
    double sum;
    for (int i = R.size() - 1; i >= 0; i--) {
        sum = 0;
        for (int j = i + 1; j < R.size(); j++) {
            sum += R[i][j] * res[j];
        }
        res[i] = (b[i] - sum) / R[i][i];
    }
    return res;
}


vector<double> InverseIterationMethod(double eigenValue, const vector<vector<double>>& A) {
    int n = A.size();
    vector<double> x_prev = { 0, 0.71, -0.61, 0.35 };
    //vector<double> x_prev(n);
    vector<vector<double>> A_copy = A;
    //x_prev[0] = 1;
    vector<double> x_next;
    vector<double> x_diff(n);
    double norm;
    int iterations = 0;

    for (int i = 0; i < n; ++i) {
        A_copy[i][i] -= eigenValue;
    }
    int prod;
    pair<vector<vector<double>>, vector<vector<double>>> QR_decomp = QRDecomposition(A_copy, prod);
    vector<vector<double>> Q = QR_decomp.first;
    vector<vector<double>> R = QR_decomp.second;
    do {
        iterations++;
        if (iterations > MAX_ITERATIONS) {
            return InverseIterationMethod(eigenValue - DELTA, A);
        }
        x_next = QR_final_step(transpose(Q), R, x_prev);
        norm = NormSquareVec(x_next);
        for (int i = 0; i < n; ++i) {
            x_next[i] /= norm;
            x_diff[i] = x_next[i] - x_prev[i];
        }
        x_prev = x_next;
    } while (NormSquareVec(matrix_prod_vec(A_copy, x_diff)) > epsilon);

    return x_next;
}


vector<vector<double>> GetEigenVectors(const vector<double>& eigen_values, vector<vector<double>> A) {
    vector<vector<double>> vectors(A.size());

    for (int i = 0; i < A.size(); ++i) {
        vectors[i] = InverseIterationMethod(eigen_values[i], A);
    }
    return vectors;
}


void EigenVectorsCheck(const vector<vector<double>>& vectors, const vector<double>& values, const vector<vector<double>>& A) {
    double discrepancy;
    vector<vector<double>> A_copy = A;
    cout << endl;
    cout << "Check eigen vectors" << endl;
    for (int i = 0; i < values.size(); ++i) {
        cout << "eigen value = " << values[i] << endl;
        cout << "vector = ";
        DisplayVector(vectors[i]);
        for (int j = 0; j < A.size(); ++j) A_copy[j][j] = A[j][j] - values[i];
        cout << "discrepancy = " << eigen_check(A_copy, vectors[i]) << endl;
        cout << endl;
    }
}


double ScalarProduct(const vector<double>& a, const vector<double>& b) {
    double res = 0;
    for (int i = 0; i < a.size(); ++i) {
        res += a[i] * b[i];
    }
    return res;
}

pair<double, vector<double>> ModifiedInverseIterationMethod(const vector<vector<double>>& A, const vector<double>& x0) {
    double lambda_next = 0, lambda_prev;
    vector<double> x_next, x_prev, x_diff(A.size());
    x_prev = x0;
    int n = A.size();
    double norm = NormSquareVec(x0);
    for (int i = 0; i < n; ++i) {
        x_prev[i] /= norm;
    }
    vector<vector<double>> A_copy = A;
    vector<vector<double>> modifiedA = A;
    int prod;
    do
    {
        A_copy = A;
        lambda_prev = lambda_next;
        lambda_next = ScalarProduct(matrix_prod_vec(A, x_prev), x_prev);
        for (int i = 0; i < n; ++i) {
            modifiedA[i][i] = A[i][i] - lambda_next;
        }
        pair<vector<vector<double>>, vector<vector<double>>> QR_decomp = QRDecomposition(modifiedA, prod);
        vector<vector<double>> Q = QR_decomp.first;
        vector<vector<double>> R = QR_decomp.second;

        x_next = QR_final_step(transpose(Q), R, x_prev);
        norm = NormSquareVec(x_next);
        for (int i = 0; i < n; ++i) {
            x_next[i] /= norm;
            A_copy[i][i] -= lambda_next;
        }
        x_prev = x_next;
    //} while (abs(lambda_next - lambda_prev) > epsilon);
    } while (NormSquareVec(matrix_prod_vec(A_copy, x_next)) > epsilon);

    cout << endl;
    cout << "ModifiedInverseIterationMethod for approximation: ";
    DisplayVector(x0);
    cout << "Eigen value = " << lambda_next << endl;
    cout << "Eigen vector: ";
    DisplayVector(x_next);
    cout << endl;
    pair<double, vector<double>> res(lambda_next, x_next);
    return res;
}






int main() {
    vector<vector<double>> Matrix;
    DataRead(Matrix, "EIGEN3.txt");
    int prodVanilaQr = 0, prodQrWithHesenberg = 0, prodQrWithShift = 0, prodQrWithHesenbergAndShift = 0;
    //DisplayVector(InverseIterationMethod(2.98, Matrix));
    //vector<double> x0 = { 0, 0.71, -0.61, 0.35 };
    //ModifiedInverseIterationMethod(Matrix, x0);
    cout << "Vanila QR decompositon" << endl;
    vector<double> eigens_v = EigenValuesQR(Matrix, prodVanilaQr);
    DisplayVector(eigens_v);
    //HesenbergDecomposition(Matrix, prodVanilaQr);
    //cout << endl;
    //cout << "QR with Hesenberg" << endl;
    //DisplayVector(EigenValuesHesenberg(Matrix, prodQrWithHesenberg));

    //cout << endl;
    //cout << "QR with Shift" << endl;
    //DisplayVector(EigenValuesShift(Matrix, prodQrWithShift));

    //cout << endl;
    //cout << "QR with Hesenberg And Shift" << endl;
    //DisplayVector(EigenValuesHesenbergAndShift(Matrix, prodQrWithHesenbergAndShift));

    //vector<vector<double>> approx_vecs = { {-1, 0.05, -0.3, -0.45},{0.1, 0.65, -0.67, 0.2}, {0.1, 0.00, -0.15, -0.1},
    //    {-0.1, 0.5, 0.25, -0.1} };
    //for (int i = 0; i < approx_vecs.size(); ++i) ModifiedInverseIterationMethod(Matrix, approx_vecs[i]);


    //cout << "Matrix A: " << endl;
    //vector<vector<double>> eigen_vectors = GetEigenVectors(eigens_v, Matrix);
    //cout << endl;
    //EigenVectorsCheck(eigen_vectors, eigens_v, Matrix);

    //cout << endl;
    //cout << "Transposed matrix A:" << endl;
    //eigen_vectors = GetEigenVectors(eigens_v, transpose(Matrix));
    //cout << endl;
    //EigenVectorsCheck(eigen_vectors, eigens_v, Matrix);
}
