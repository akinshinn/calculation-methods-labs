#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
using namespace std;
const double epsilon = 1e-4;
const double epsilonzero = 1e-4;
const double omega = 0.5;
void DataRead(int& SystemNumber, vector<int>& size, vector<vector<vector<double>>>& Matrix, string filename) {
    ifstream file(filename);
    string line;
    getline(file, line);
    stringstream str(line);
    int intSize = 0;
    vector<vector<double>> mas{};
    str >> SystemNumber;
    for (int k = 0; k < SystemNumber; k++) {
        getline(file, line);
        stringstream str(line);
        str >> intSize;
        mas.resize(intSize, vector<double>(intSize + 1, 0));
        for (auto& row : mas) {
            row.resize(intSize + 1, 0);
        }

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

void Display2DMatrix(vector<vector<double>>& Matrix) {
    for (auto& row : Matrix) {
        for (auto& el : row) {
            cout << el << " ";
        }
        cout << endl;
    }
}

double EuclideNorm(vector<double>& vec1, vector<double>& vec2) {
    double suma = 0;
    for (int i = 0; i < vec1.size(); i++) {
        suma += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return suma;
}

bool StopCriteriaOne(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol) {

    return EuclideNorm(CurIterSol, NextIterSol) < epsilon;
}

bool StopCriteriaSecond(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol) {
    vector<double> nul{};
    nul.resize(size, 0);
    return EuclideNorm(CurIterSol, NextIterSol) / (EuclideNorm(CurIterSol, nul) + epsilonzero) < epsilon;
}

bool StopCriteriaThird(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol) {
    vector<double> f{}, eval{};
    f.resize(size, 0);
    eval.resize(size, 0);
    for (int i = 0; i < size; i++) {
        f[i] = matrix[i][size];
    }
    double el = 0;
    for (int i = 0; i < size; i++) {
        el = 0;
        for (int j = 0; j < size; j++) {
            el += matrix[i][j] * CurIterSol[j];
        }
        eval[i] = el;
    }
    return EuclideNorm(f, eval) < epsilon;
}


pair<vector<double>, int> JacobyMethod(int size, vector<vector<double>>& matrix, bool(&StopCriteria)(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol)) {
    vector<double> CurIterSol{}, NextIterSol{};
    double suma = 0;
    int iteration = 0;
    CurIterSol.resize(size, 0);
    NextIterSol.resize(size, 0);
    do {
        iteration++;
        CurIterSol = NextIterSol;
        for (int i = 0; i < size; i++) {
            suma = 0;
            for (int j = 0; j < i; j++) {
                suma += matrix[i][j] * CurIterSol[j];
            }
            for (int j = i + 1; j < size; j++) {
                suma += matrix[i][j] * CurIterSol[j];
            }
            NextIterSol[i] = (matrix[i][size] - suma) / matrix[i][i];
        }
    } while (!StopCriteria(size, matrix, CurIterSol, NextIterSol));
    // StopCriteria - true если надо остановиться
    return { NextIterSol , iteration };
}

void WriteJacobyAnswer(int SystemNumber, vector<int>& size, vector<vector<vector<double>>>& Matrix, string filename) {
    ofstream out;
    pair<vector<double>, int> res;
    out.open(filename);
    out << "JacobyMethod" << endl;
    out << "epsilon = " << epsilon << endl;
    out << endl;

    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = JacobyMethod(size[i], Matrix[i], StopCriteriaOne);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }
    out << endl;
    out << "StopCriteria ||x^(k+1)-x^k|| / (||x^k||+epsilonzero) < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = JacobyMethod(size[i], Matrix[i], StopCriteriaSecond);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out << endl;
    out << "StopCriteria ||Ax^k-f|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = JacobyMethod(size[i], Matrix[i], StopCriteriaThird);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out.close();
}


pair<vector<double>, int> RelaxationMethod(int size, double w, vector<vector<double>>& matrix, bool(&StopCriteria)(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol)) {
    vector<double> CurIterSol{}, NextIterSol{};
    double suma = 0;
    int iteration = 0;
    CurIterSol.resize(size, 0);
    NextIterSol.resize(size, 0);
    do {
        iteration++;
        CurIterSol = NextIterSol;
        for (int i = 0; i < size; i++) {
            suma = 0;
            for (int j = 0; j < i; j++) {
                suma += matrix[i][j] * NextIterSol[j];
            }
            for (int j = i + 1; j < size; j++) {
                suma += matrix[i][j] * CurIterSol[j];
            }
            NextIterSol[i] = (1 - w) * CurIterSol[i] + w * (matrix[i][size] - suma) / matrix[i][i];
        }
    } while (!StopCriteria(size, matrix, CurIterSol, NextIterSol));
    // StopCriteria - true если надо остановиться
    return { NextIterSol , iteration };
}


void WriteRelaxationAnswer(int SystemNumber, vector<int>& size, vector<vector<vector<double>>>& Matrix, string filename) {
    ofstream out;
    pair<vector<double>, int> res;
    out.open(filename);
    out << "RelaxationMethod" << endl;
    out << "epsilon = " << epsilon << endl;
    out << endl;

    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = RelaxationMethod(size[i], omega, Matrix[i], StopCriteriaOne);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }
    out << endl;
    out << "StopCriteria ||x^(k+1)-x^k|| / (||x^k||+epsilonzero) < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = RelaxationMethod(size[i], omega, Matrix[i], StopCriteriaSecond);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out << endl;
    out << "StopCriteria ||Ax^k-f|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = RelaxationMethod(size[i], omega, Matrix[i], StopCriteriaThird);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out.close();
}


void main()
{
    int SystemNumber;
    vector<int> size{};
    vector<vector<vector<double>>> Matrix{};

    DataRead(SystemNumber, size, Matrix, "System.txt");

    WriteJacobyAnswer(SystemNumber, size, Matrix, "JacobyAnswer.txt");

    WriteRelaxationAnswer(SystemNumber, size, Matrix, "RelaxationAnswer.txt");
}