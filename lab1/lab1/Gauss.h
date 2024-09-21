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
vector<float> readf_vector_float(string file);




vector<vector<double>> readf_matrix_double(string file_matrix);

vector<vector<float>> readf_matrix_float(string file_matrix);

template <typename T>
bool check_matrix(const vector<vector<T>>& triang_matrix) {
    T num = 1;
    for (int i = 0; i < triang_matrix.size(); ++i) {
        num *= triang_matrix[i][i];
    }
    return (abs(num) > 1e-15);
}

template <typename T>
void triangularize(vector<vector<T>>& matrix) {
    int n = matrix.size();

    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            T temp = matrix[i][j];
            for (int k = 0; k < n + 1; ++k) {
                matrix[i][k] -= matrix[j][k] * temp / matrix[j][j];
            }
        }
    }

}

template <typename T>
void print_matrix(const vector<vector<T>>& matrix) {
    int n = matrix.size();
    cout << n << endl;
    if (matrix.size() < matrix[0].size()) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n + 1; ++j) {
                cout << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }
    else {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }

    cout << endl;
}



vector<vector<double>> read_matrix();





template <typename T>
int get_main_element_row(const vector<vector<T>>& matrix, int var_row) {
    T max = abs(matrix[var_row][var_row]);
    int max_row = var_row;
    int n = matrix.size();
    for (int i = var_row + 1; i < n; ++i) {
        if (matrix[i][var_row] > max) {
            max = abs(matrix[i][var_row]);
            max_row = i;
        }
    }
    return max_row;
}


template <typename T>
void permutate_rows(vector<vector<T>>& matrix) {
    int n = matrix.size();
    int main_elem_row;
    vector<T> temp_vec = matrix[0];
    for (int i = 0; i < n; ++i) {
        main_elem_row = get_main_element_row(matrix, i);
        temp_vec = matrix[main_elem_row];
        matrix[main_elem_row] = matrix[i];
        matrix[i] = temp_vec;
    }
}

template <typename T>
vector<T> Gauss_method(vector<vector<T>>& matrix) {
    //permutate_rows(matrix);
    print_matrix(matrix);
    triangularize(matrix);
    print_matrix(matrix);
    int n = matrix.size();
    vector<T> x(n);
    if (!check_matrix(matrix)) {
        cout << "det A = 0" << endl;
        return x;
    }

    T sum;
    for (int i = n - 1; i > -1; --i) {
        sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += matrix[i][j] * x[j];
        }
        x[i] = (matrix[i][n] - sum) / matrix[i][i];
        cout << x[i] << endl;
    }
    return x;
}


template <typename T>
vector<T> check_ans(const vector<vector<T>>& matrix, const vector<T>& x) {
    int n = matrix.size();
    vector<T> diff(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            diff[i] += matrix[i][j] * x[j];
        }
        diff[i] -= matrix[i][n];
    }
    return diff;
}




template <typename T>
vector<vector<T>> matrix_b_concat(const vector<vector<T>>& matrix, const vector<T>& b) {
    int n = matrix.size();
    vector<vector<T>> res;

    for (int i = 0; i < n; ++i) {
        res.emplace_back(vector<T>(n + 1));
        for (int j = 0; j < n; ++j) {
            res[i][j] = matrix[i][j];
        }
        res[i][n] = b[i];
    }
    return res;
}



template <typename T>
void write_file(string file,const vector<T>& x){
    ofstream out;
    out.open(file);
    //T discrepancy = Discrepancy(x,)
    if (out.is_open()) {
        out << "solution" << endl;
        for (int i = 0; i < x.size(); ++i) {
            out << x[i] << " ";
        }
        out.close();
    }
}


const float fsqrt = 0.5;
const float fsquare = 2;
const double dsqrt = 0.5;
const double dsquare = 2;
const float outrageFloat = 0.1;
const double outrageDouble = 0.1;
const double flagOfNoSolution = -10;

pair<vector<float>, int> QRFloat(int size, vector<vector<float>>& Matrix);
pair<vector<float>, int> QRFloatDecomposition(int size, vector<vector<float>>& Matrix, vector<vector<float>>& TMatrix);
pair<vector<double>, int> QRDouble(int size, vector<vector<double>>& Matrix);
pair<vector<double>, int> QRDoubleDecomposition(int size, vector<vector<double>>& Matrix, vector<vector<double>>& TMatrix);
vector<vector<float>> InverseFloat(vector<vector<float>> Matrix);
void FloatTest(string filename);
void DoubleTest(string filename);
template<typename T>
T NormOneMatrix(vector<vector<T>>& mas) {
    T suma = 0;
    T MaxSuma = -numeric_limits<T>::max();
    for (int j = 0; j < mas[0].size(); j++) {
        suma = 0;
        for (int i = 0; i < mas.size(); i++) {
            suma += abs(mas[i][j]);
        }
        MaxSuma = max(MaxSuma, suma);

    }
    return MaxSuma;
}

template<typename T>
T NormInfMatrix(vector<vector<T>>& mas) {
    T suma = 0;
    T MaxSuma = -numeric_limits<T>::max();
    for (int i = 0; i < mas.size(); i++) {
        suma = 0;
        for (int j = 0; j < mas[0].size(); j++) {
            suma += abs(mas[i][j]);
        }
        MaxSuma = max(MaxSuma, suma);
    }
    return MaxSuma;
}


template<typename T>
T NormSquareVec(vector<T>& vec) {
    T suma = 0;
    for (int i = 0; i < vec.size(); i++) {
        suma += vec[i] * vec[i];
    }
    return sqrt(suma);
}


template<typename T>
T NormOneVec(vector<T>& vec) {
    T suma = 0;
    for (int i = 0; i < vec.size(); i++) {
        suma += abs(vec[i]);
    }
    return suma;
}

template<typename T>
T NormInfVec(vector<T>& vec) {
    T Maxel = -numeric_limits<T>::max();
    for (int i = 0; i < vec.size(); i++) {
        Maxel = max(Maxel, abs(vec[i]));
    }
    return Maxel;
}


template<typename T>
vector<vector<T>> MultiplyMatrix(vector<vector<T>>& mas1, vector<vector<T>>& mas) {
    vector<vector<T>> res{};
    int size = mas1.size();
    res.resize(size, vector<T>(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                res[i][j] += mas1[i][k] * mas[k][j];
            }
        }
    }
    return res;
}


template<typename T>
void DisplayAnswer(vector<vector<T>>& ans) {
    int example = 1;
    int i = 1;
    for (auto& sol : ans) {
        i = 1;
        cout << "Exmaple " << example++ << endl;
        if (sol.second == 1) {
            for (auto& el : sol.first) {
                cout << "x" << i++ << "=" << el << endl;
            }
        }
        else {
            cout << "The system has infinitely many solutions";
        }
    }

}


template<typename T>
void WriteAnswer(vector<pair<vector<T>, int>>& ans, string filename) {
    ofstream out;
    out.open(filename);
    int example = 1;
    int i = 1;
    for (auto& sol : ans) {
        i = 1;
        out << "Exmaple " << example++ << endl;
        if (sol.second == 1) {
            for (auto& el : sol.first) {
                out << "x" << i++ << "=" << el << endl;
            }
        }
        else {
            out << "Матрица системы вырожденная, система несовместна" << endl;
        }
    }
    out.close();
}


template<typename T>
void DisplayMatrix(vector<vector<T>>& Matrix) {
    int size = Matrix.size();
    for (int i = 0; i < size; i++) {
        for (int j = 0; j <= size; j++) {
            cout << Matrix[i][j] << " ";
        }
        cout << endl;
    }
}

template<typename T>
void DisplayVector(vector<T>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        cout << "x" << i + 1 << "=" << vec[i] << endl;
    }
}


template<typename T>
void DataRead(int& SystemNumber, vector<int>& size, vector<vector<vector<T>>>& Matrix, string filename) {
    ifstream file(filename);
    string line;
    getline(file, line);
    stringstream str(line);
    int intSize = 0;
    vector<vector<T>> mas{};
    str >> SystemNumber;
    for (int k = 0; k < SystemNumber; k++) {
        getline(file, line);
        stringstream str(line);
        str >> intSize;
        mas.resize(intSize, vector<T>(intSize + 1));
        for (int Curstr = 0; Curstr < intSize; Curstr++) {
            getline(file, line);
            stringstream str(line);
            for (int j = 0; j <= intSize; j++) {
                str >> mas[Curstr][j];
            }
        }
        size.push_back(intSize);
        Matrix.push_back(mas);
    }

    file.close();
}

template<typename T>
T Discrepancy(vector<T>& vec1, vector<T>& vec2) {
    T ans = 0;
    T suma = 0;
    if (vec1.size() == vec2.size()) {
        for (int i = 0; i < vec1.size(); i++) {
            suma = suma + (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
        }
        ans = sqrt(suma);
    }
    else {
        cout << "Разные размерности у векторов, для которых ищем невязку";
    }
    return ans;
}



template<typename T>
vector<T> SearchRightParts(vector<T>& solution, vector<vector<T>>& SystemMatrix) {
    T suma = 0;
    vector<T> CalculatedRightParts{};
    if (solution.size() == SystemMatrix.size()) {
        for (int i = 0; i < solution.size(); i++) {
            suma = 0;
            for (int j = 0; j < solution.size(); j++) {
                suma += SystemMatrix[i][j] * solution[j];
            }
            CalculatedRightParts.push_back(suma);
        }
    }
    else {
        cout << "Проблема размерности в SearchRightParts";
    }
    return CalculatedRightParts;
}



template<typename T>
void WriteDiscrepancy(vector<T>& ans, string filename) {
    ofstream out;
    out.open(filename);
    out << "Невязка" << endl;
    int example = 1;
    int i = 1;
    for (auto& el : ans) {
        if (el == flagOfNoSolution) {
            out << "Example " << example++ << ": Матрицы системы вырождена" << endl;
        }
        else {
            out << "Example " << example++ << ": " << el << endl;
        }
    }
    out.close();
}


template<typename T>
void WriteConditionality(vector<T>& masOne, vector<T>& masInfty, string filename) {
    ofstream out;
    out.open(filename);
    for (int i = 0; i < masOne.size(); i++) {
        if (masOne[i] != flagOfNoSolution) {
            out << "Example" << i + 1 << endl << " NormOne: " << masOne[i] << ", NormInfty: " << masInfty[i] << endl;
        }
        else {
            out << "Example" << i + 1 << endl << " Матрица несовместна" << endl;
        }
    }

    out.close();
}


template<typename T>
void WriteCheckPoint(vector<vector<vector<T>>>& mas, string filename) {
    ofstream out;
    out.open(filename);
    int size = 0, SystemNumbers = mas.size();
    for (int k = 0; k < SystemNumbers; k++) {
        size = mas[k].size();
        out << "Example: " << k + 1 << endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                out << mas[k][i][j] << " ";
            }
            out << endl;
        }
    }
    out.close();
}


template<typename T>
void WriteAssessmentConditionality(vector<T>& vecOne, vector<T>& vecInfty, string filename) {
    ofstream out;
    out.open(filename);
    out << "AssessmentConditionality" << endl;
    for (int i = 0; i < vecOne.size(); i++) {
        if (vecOne[i] != flagOfNoSolution) {
            out << "Example: " << i + 1 << endl << "NormOne: " << vecOne[i] << ", NormInfty: " << vecInfty[i] << endl;

        }
        else {
            out << "Example: " << i + 1 << endl << "Матрица не совместна" << endl;
        }
    }
    out.close();
}

template<typename T>
void WriteQRMatrix(vector<vector<vector<T>>>& mas1, vector<vector<vector<T>>>& mas2, string filename) {
    ofstream out;
    out.open(filename);
    out << "QR decomposition" << endl;
    for (int i = 0; i < mas1.size(); i++) {
        if (mas1[i].size() != 0) {
            out << "Example: " << i + 1 << endl;
            out << "Q Matrix" << endl;
            for (int j = 0; j < mas2[i].size(); j++) {
                for (int u = 0; u < mas2[i].size(); u++) {
                    out << mas2[i][u][j] << " ";
                }
                out << endl;
            }
            out << "R Matrix" << endl;
            for (auto& row : mas1[i]) {
                for (int u = 0; u < row.size() - 1; u++) {
                    out << row[u] << " ";
                }
                out << endl;
            }

        }
        else {
            out << "Example: " << i + 1 << endl << "Матрица не совместна" << endl;
        }
    }
    out.close();
}


template<typename T>
T AssessmentConditionality(vector<T>& trueSol, vector<T>& OutragedSol, vector<T>& trueRigthPart, vector<T>& OutragedRigthPart, T(&norm)(vector<T>&)) {
    vector<T> deltaSol{}, deltaRigthPart{};
    int size = trueSol.size();
    for (int i = 0; i < size; i++) {
        deltaSol.push_back(trueSol[i] - OutragedSol[i]);
        deltaRigthPart.push_back(trueRigthPart[i] - OutragedRigthPart[i]);
    }
    T normDeltaSol = norm(deltaSol), normSol = norm(trueSol), normDeltaRightPart = norm(deltaRigthPart), normRigthPart = norm(trueRigthPart);
    return normDeltaSol * normRigthPart / (normSol * normDeltaRightPart);
}


template<typename T>
T SearchConditionality(vector<vector<T>> mas1, vector<vector<T>> mas2, T(&norm)(vector<vector<T>>&)) {
    if (mas1.size() != mas1[0].size()) {
        for (int i = 0; i < mas1.size(); i++) {
            mas1[i].pop_back();
        }
    }
    if (mas2.size() != mas2[0].size()) {
        for (int i = 0; i < mas2.size(); i++) {
            mas2[i].pop_back();
        }
    }
    return norm(mas1) * norm(mas2);
}


//void main2()
//{
//    setlocale(LC_ALL, "Russian");
//    DoubleTest("LinearSystem.txt");
//
//}