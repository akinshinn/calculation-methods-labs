#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
using namespace std;
const double epsilon = 1e-5;
const double epsilonzero = 1e-4;
const double omega = 0.5;
const double tau = 1e-2;
const int maxIterations = 100;

void Display2DMatrix(vector<vector<double>>& Matrix);
void DisplayVector(const vector<double> x);
vector<vector<double>> matrix_prod(const vector<vector<double>>& A, const vector<vector<double>>& B);


double NormSquareVec(const vector<double>& vec) {
    double suma = 0;
    for (int i = 0; i < vec.size(); i++) {
        suma += vec[i] * vec[i];
    }
    return sqrt(suma);
}


double NormOneVec(const vector<double>& vec) {
    double suma = 0;
    for (int i = 0; i < vec.size(); i++) {
        suma += abs(vec[i]);
    }
    return suma;
}

double NormInfVec(const vector<double>& vec) {
    double Maxel = -numeric_limits<double>::max();
    for (int i = 0; i < vec.size(); i++) {
        Maxel = max(Maxel, abs(vec[i]));
    }
    return Maxel;
}


double NormInfMatrix(vector<vector<double>>& mas) {
    double suma = 0;
    double MaxSuma = -1;
    for (int i = 0; i < mas.size(); i++) {
        suma = 0;
        for (int j = 0; j < mas.size(); j++) {
            suma += abs(mas[i][j]);
        }
        MaxSuma = max(MaxSuma, suma);
    }
    return MaxSuma;
}


double NormOneMatrix(vector<vector<double>>& mas) {
    double suma = 0;
    double MaxSuma = -numeric_limits<double>::max();
    for (int j = 0; j < mas.size(); j++) {
        suma = 0;
        for (int i = 0; i < mas.size(); i++) {
            suma += abs(mas[i][j]);
        }
        MaxSuma = max(MaxSuma, suma);

    }
    return MaxSuma;
}


vector<vector<double>> get_reverse_LD_matrix( vector<vector<double>> A) {
    vector<vector<double>> LD_rev,inverse;
    int size = A.size();
    vector<double> x(size);

    for (int i = 0; i < size; i++) {
        A[i][size] = 0;
    }
    for (int i = 0; i < size; i++) {
        A[i][size] = 1;
        if (i != 0) A[i - 1][size] = 0;
        

        for (int j = 0; j < size; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += A[j][k] * x[k];
            }
            x[j] = (A[j][size] - sum) / A[j][j];
        }
        
        inverse.push_back(x);

    }
    LD_rev.resize(size, vector<double>(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            LD_rev[i][j] = inverse[j][i];
        }
    }
    return LD_rev;
}


// Первое значение - норма C, второе - норма C_L, третье - норма C_U
vector<double> getNormC_SIM(const vector<vector<double>>& A, double tau, double(&norm)(vector<vector<double>>& matrix) = NormInfMatrix) {
    int n = A.size();
    vector<vector<double>> C;
    C.resize(n, vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                C[i][j] = -tau * A[i][j];

            }
            else {
                C[i][i] = -tau * A[i][i] + 1;
            }
        }
    }

    vector<vector<double>> C_L = C, C_U = C;
    for (int i = 0; i < n; ++i) {
        C_L[i][i] = 0;
        C_U[i][i] = 0;
        for (int j = 0; j < n; ++j) {
            if (j < i) {
                C_U[i][j] = 0;
            }
            else if (j > i) {
                C_L[i][j] = 0;
            }

        }
    }
    vector<double> res= {norm(C), norm(C_L), norm(C_U)};
    return res;
}

 //Первое значение - норма C, второе - норма C_L, третье - норма C_U
vector<double> getNormC_Seidel(const vector<vector<double>>& A, double(&norm)(vector<vector<double>>& matrix) = NormOneMatrix){
    vector<vector<double>> C;
    int size = A.size();
    vector<vector<double>> LD_reverse = get_reverse_LD_matrix(A);
    C.resize(size, vector<double>(size));
    vector<vector<double>> prod = matrix_prod(LD_reverse, A);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != j) {
                C[i][j] = -prod[i][j];
            }
            else {
                C[i][i] = 1 - prod[i][i];
            }
        }
    }

    vector<vector<double>> C_L = C, C_U = C;
    for (int i = 0; i < size; ++i) {
        C_L[i][i] = 0;
        C_U[i][i] = 0;
        for (int j = 0; j < size; ++j) {
            if (j < i) {
                C_U[i][j] = 0;
            }
            else if (j > i) {
                C_L[i][j] = 0;
            }

        }
    }
    vector<double> res = { norm(C), norm(C_L), norm(C_U) };
    return res;
}


vector<double> vec_diff(const vector<double>& v1, const vector<double>& v2) {
    vector<double> res(v1.size());
    for (int i = 0; i < v1.size(); ++i) {
        res[i] = v1[i] - v2[i];
    }
    return res;
}

double tau_estimate(const vector<vector<double>>& matrix) {
    double maxx = -1;
    bool flag = 1;
    double sum;
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        sum = 0;
        for (int j = 0; i < i; ++j) {
            sum += abs(matrix[i][j]);
        }
        for (int j = i + 1; j < n; ++j) {
            sum += abs(matrix[i][j]);
        }
        maxx = max(sum + abs(matrix[i][i]), maxx);
        if (sum > abs(matrix[i][i])) {
            flag = 0;
            break;
        }
    }
    if (flag) {
        return maxx;
    }
    return -1;
}


vector<double> matrix_prod_vec(const vector<vector<double>>& A, const vector<double> b) {
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

vector<vector<double>> matrix_prod(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int n = A.size();
    vector<vector<double>> res;
    res.resize(n, vector<double>(n));


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += A[i][k] * B[k][j];
            }
            res[i][j] = sum;

        }

    }
    return res;
}


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


void DisplayVector(const vector<double> x) {
    cout << "Size: " << x.size() << endl;
    for (int i = 0; i < x.size(); ++i) {
        cout << x[i] << " ";
    }
    cout << endl;
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
    return sqrt(suma);
}
//bool(&StopCriteria)(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol)
bool StopCriteriaOne(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol) {
    //cout << EuclideNorm(CurIterSol, NextIterSol) << endl;
    return EuclideNorm(CurIterSol, NextIterSol) < epsilon;
}

bool StopCriteriaSecond(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol) {
    vector<double> nul{};
    nul.resize(size, 0);
    //cout << EuclideNorm(CurIterSol, NextIterSol) / (EuclideNorm(CurIterSol, nul) + epsilonzero) << endl;    
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
    //cout << EuclideNorm(f, eval) << endl;
    return EuclideNorm(f, eval) < epsilon;
}

bool StopCriteriaThirdTriangle(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol) {
    vector<double> f{}, eval{};
    f.resize(size, 0);
    eval.resize(size, 0);
    for (int i = 0; i < size; i++) {
        f[i] = matrix[i][3];
    }
    double el = 0;
    double a, b, c;
    for (int i = 0; i < size; i++) {
        a = matrix[i][0];
        b = matrix[i][1];
        c = matrix[i][2];

        el = a * CurIterSol[max(i - 1,0)] + b * CurIterSol[i] + c * CurIterSol[min(i + 1, size-1)];
        eval[i] = el;
    }
    //cout << EuclideNorm(f, eval) << endl;
    return EuclideNorm(f, eval) < epsilon;
}


bool StopCriteriaFour(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol)
{
    vector<double> x0(size);
    vector<double> b(size);
    for (int i = 0; i < size; ++i) {
        b[i] = matrix[i][size];
    }
    vector<double> temp = matrix_prod_vec(matrix, x0);
    double r0 = EuclideNorm(b, temp);
    vector<double> temp2 = matrix_prod_vec(matrix, NextIterSol);
    double rk = EuclideNorm(b, temp2);
    return (rk / r0) < epsilon;
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
        if (iteration > 200) break;
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

void GenerateThreeDiagonal(int n, string filename) {
    ofstream out;
    pair<vector<double>, int> res;
    out.open(filename);
    out << n << endl;
    out << "0 4 1 6" << endl;
    for (int i = 2; i < n; i++) {
        out << "1 4 1 " << 10 - 2 * (i % 2) << endl;
    }
    out << "1 4 0 " << 9 - 3 * (n % 2) << endl;
    out.close();
}

void DataReadTriangleMatrix(int& SizeTriangleMatrix, vector<vector<double>>& TriangleMatrix, string filename) {
    ifstream file(filename);
    string line;
    getline(file, line);
    stringstream str(line);
    vector<double> mas{};
    mas.resize(4, 0);
    str >> SizeTriangleMatrix;
    for (int k = 0; k < SizeTriangleMatrix; k++) {
        getline(file, line);
        stringstream str(line);
        for (int j = 0; j < 4; j++) {
            str >> mas[j];
        }
        TriangleMatrix.push_back(mas);
    }

    file.close();
}


pair<vector<double>, int> RelaxationMethodTriangle(int size, double w, vector<vector<double>>& matrix, bool(&StopCriteria)(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol)) {
    vector<double> CurIterSol{}, NextIterSol{};
    double suma = 0, a, b, c, d;
    int iteration = 0;
    CurIterSol.resize(size, 0);
    NextIterSol.resize(size, 0);
    do {
        iteration++;
        CurIterSol = NextIterSol;
        for (int i = 0; i < size; i++) {
            a = matrix[i][0];
            b = matrix[i][1];
            c = matrix[i][2];
            d = matrix[i][3];
            NextIterSol[i] = (1 - w) * CurIterSol[i] + w * (d - a*NextIterSol[max(i-1,0)] - c * CurIterSol[min(i+1, size-1)]) / b;
        }
    } while (!StopCriteria(size, matrix, CurIterSol, NextIterSol));
    // StopCriteria - true если надо остановиться
    return { NextIterSol , iteration };
}


void WriteRelaxationTriangleAnswer(int size, vector<vector<double>>& Matrix, string filename) {
    ofstream out;
    pair<vector<double>, int> res;
    out.open(filename);
    out << "RelaxationMethod" << endl;
    out << "epsilon = " << epsilon << endl;
    out << endl;

    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    res = RelaxationMethodTriangle(size, omega, Matrix, StopCriteriaOne);
    out << "count of iteration = " << res.second << endl;
    for (int j = 0; j < size; j++) {
        out << "x" << j + 1 << "=" << res.first[j] << " ";
    }
    out << endl;

    out << "StopCriteria ||x^(k+1)-x^k|| / (||x^k||+epsilonzero) < epsilon" << endl;
    res = RelaxationMethodTriangle(size, omega, Matrix, StopCriteriaSecond);
    out << "count of iteration = " << res.second << endl;
    for (int j = 0; j < size; j++) {
        out << "x" << j + 1 << "=" << res.first[j] << " ";
    }
    out << endl;


    out << "StopCriteria ||Ax^k-f|| < epsilon" << endl;
    res = RelaxationMethodTriangle(size, omega, Matrix, StopCriteriaThirdTriangle);
    out << "count of iteration = " << res.second << endl;
    for (int j = 0; j < size; j++) {
        out << "x" << j + 1 << "=" << res.first[j] << " ";
    }
    out << endl;

    out.close();
}


void OmegaVsIteration(int& SystemNumber, vector<int>& size, vector<vector<vector<double>>>& matrix, string filename, string Wolfram) {
    ofstream out, wolf;
    out.open(filename);
    wolf.open(Wolfram);
    pair<vector<double>, int> res;
    for (double w = 0.1; w < 2; w += 0.1) {
        wolf << w << " ";
    }
    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        out << "System :" << i + 1 << endl;
        wolf << endl;
        for (double w = 0.1; w < 2; w += 0.1) {
            res = RelaxationMethod(size[i], w, matrix[i], StopCriteriaOne);
            out << "w = " << w << ", iteration = " << res.second << endl;
            wolf << res.second << " ";
        }
    }
    out << endl;
    out << "StopCriteria ||x^(k+1)-x^k|| / (||x^k||+epsilonzero) < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        out << "System :" << i + 1 << endl;
        wolf << endl;
        for (double w = 0.1; w < 2; w += 0.1) {
            res = RelaxationMethod(size[i], w, matrix[i], StopCriteriaSecond);
            out << "w = " << w << ", iteration = " << res.second << endl;
            wolf << res.second << " ";
        }
    }
    out << endl;
    out << "StopCriteria ||Ax^k-f|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        out << "System :" << i + 1 << endl;
        
        wolf << endl;
        for (double w = 0.1; w < 2; w += 0.1) {
            res = RelaxationMethod(size[i], w, matrix[i], StopCriteriaThird);
            out << "w = " << w << ", iteration = " << res.second << endl;
            wolf << res.second << " ";
        }
    }
    out << endl;
    out.close();
    wolf.close();
}


pair<vector<double>, int> SimpleIterationMethod(int size, vector<vector<double>>& matrix, bool(&StopCriteria)(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol),
    pair<bool, double> tau_var = {0, 0}) {
    vector<double> nextSol(size), curSol(size);

    vector<double> temp;
    int iteration = 0;
    //Display2DMatrix(matrix);
    double tau_cur;
    if (tau_var.first) {
        tau_cur = tau_var.second;
    }
    else {
        tau_cur = 1 / (tau_estimate(matrix));

        if (tau_cur < 0)
        {
            tau_cur = tau;
        }
    }
    do {
        iteration++;
        curSol = nextSol;
        temp = matrix_prod_vec(matrix, curSol);
        for (int i = 0; i < size; ++i) {
            nextSol[i] = curSol[i] - tau_cur * temp[i] + tau_cur * matrix[i][size];
        }
        if (iteration >= maxIterations) break;
        cout << tau_cur;
        //DisplayVector(nextSol);
        
    } while (!StopCriteria(size, matrix, curSol, nextSol));
    return { nextSol, iteration };
  }


pair<vector<double>, int> SeidelMethod(int size, vector<vector<double>>& matrix, bool(&StopCriteria)(int size, vector<vector<double>>& matrix, vector<double>& CurIterSol, vector<double>& NextIterSol)) {
    //vector<double> nextSol(size), curSol(size);

    vector<double> temp, nextSol = {0,0,0,0};
    vector<double> curSol;
    int iteration = 0;
    do {
        iteration++;
        curSol = nextSol;
        for (int j = 0; j < size; ++j) {
            double s1 = 0, s2 = 0;

            for (int i = 0; i < j; ++i) {
                s1 += matrix[j][i] * nextSol[i];
            }
            for (int i = j+1; i < size; ++i) {
                s2 += matrix[j][i] * curSol[i];
            }

            nextSol[j] = (matrix[j][size] - s1 - s2) / matrix[j][j];
        }

    } while (!StopCriteria(size, matrix, curSol, nextSol));
    return { nextSol, iteration };
}


void WriteSimpleIterationAnswer(int SystemNumber, vector<int>& size, vector<vector<vector<double>>>& Matrix, string filename) {
    ofstream out;
    pair<vector<double>, int> res;
    out.open(filename);

    out << "Simple Iteration Method" << endl;
    out << "epsilon = " << epsilon << endl;
    out << endl;
    out << endl;
    for (int i = 0; i < SystemNumber; ++i) {
        out << "Matrix " << i + 1 << endl;
        double tau_cur = 1 / (tau_estimate(Matrix[i]));
        if (tau_cur < 0) tau_cur = tau;
        vector<double> norms = getNormC_SIM(Matrix[i], tau_cur);
        out << "Norm C = " << norms[0] << endl;
        out << "Norm C_L = " << norms[1] << endl;
        out << "Norm C_U = " << norms[2] << endl;
    }
    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = SimpleIterationMethod(size[i], Matrix[i], StopCriteriaOne);

        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;

        double tau_cur = 1 / (tau_estimate(Matrix[i]));
        if (tau_cur < 0) tau_cur = tau;
        out << "tau = " << tau_cur << endl;

        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
        
    }
    out << endl;
    out << "StopCriteria ||x^(k+1)-x^k|| / (||x^k||+epsilonzero) < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = SimpleIterationMethod(size[i], Matrix[i], StopCriteriaSecond);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;

        double tau_cur = 1 / (tau_estimate(Matrix[i]));
        if (tau_cur < 0) tau_cur = tau;
        out << "tau = " << tau_cur << endl;

        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out << endl;
    out << "StopCriteria ||Ax^k-f|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = SimpleIterationMethod(size[i], Matrix[i], StopCriteriaThird);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;

        double tau_cur = 1 / (tau_estimate(Matrix[i]));
        if (tau_cur < 0) tau_cur = tau;
        out << "tau = " << tau_cur << endl;

        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out.close();
}


void WriteSeidelAnswer(int SystemNumber, vector<int>& size, vector<vector<vector<double>>>& Matrix, string filename) {
    ofstream out;
    pair<vector<double>, int> res;
    out.open(filename);
    out << "Seidel Method" << endl;
    out << "epsilon = " << epsilon << endl;
    out << endl;
    for (int i = 0; i < SystemNumber; ++i) {
        out << "Matrix " << i + 1 << endl;
        vector<double> norms = getNormC_Seidel(Matrix[i]);
        out << "Norm C = " << norms[0] << endl;
        out << "Norm C_L = " << norms[1] << endl;
        out << "Norm C_U = " << norms[2] << endl;
    }

    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    for (int i = 0; i < SystemNumber; i++) {
        res = SeidelMethod(size[i], Matrix[i], StopCriteriaOne);
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
        res = SeidelMethod(size[i], Matrix[i], StopCriteriaSecond);
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
        res = SeidelMethod(size[i], Matrix[i], StopCriteriaThird);
        out << "Example " << i + 1 << endl;
        out << "count of iteration = " << res.second << endl;
        for (int j = 0; j < size[i]; j++) {
            out << "x" << j + 1 << "=" << res.first[j] << " ";
        }
        out << endl;
    }

    out.close();
}


void TauVsIteration(vector<vector<double>>& matrix, string filename, string Wolfram, int k = 100) {
    vector<double> taus(k);
    taus[0] = 0.01;
    vector<int> iterations1(k), iterations2(k), iterations3(k);
    for (int i = 1; i < k; ++i) {
        taus[i] =taus[i - 1] + (double)1/k;
        iterations1[i] = SimpleIterationMethod(matrix.size(), matrix, StopCriteriaOne,{1,taus[i]}).second;
        iterations2[i] = SimpleIterationMethod(matrix.size(), matrix, StopCriteriaSecond, { 1,taus[i] }).second;
        iterations3[i] = SimpleIterationMethod(matrix.size(), matrix, StopCriteriaThird, { 1,taus[i] }).second;
    }

    ofstream out;
    out.open(filename);

    out << "StopCriteria ||x^(k+1)-x^k|| < epsilon" << endl;
    for (int i = 0; i < k; i++) {

        out << "count of iteration = " << iterations1[i] << ", ";

        out << "tau = " << taus[i] << endl;

        out << endl;
    }
    out << endl;
    out << "StopCriteria ||x^(k+1)-x^k|| / (||x^k||+epsilonzero) < epsilon" << endl;
    for (int i = 0; i < k; i++) {

        out << "count of iteration = " << iterations2[i] << ", ";

        out << "tau = " << taus[i] << endl;

        out << endl;
    }
    out << endl;
    out << "StopCriteria ||Ax^k-f|| < epsilon" << endl;
    for (int i = 0; i < k; i++) {

        out << "count of iteration = " << iterations3[i] << ", ";

        out << "tau = " << taus[i] << endl;

        out << endl;
    }
    out.close();
    out.open(Wolfram);
    out << k << endl;
    for (int i = 0; i < k; i++) {
        out<<  taus[i] << " " <<  iterations1[i] <<   endl;
    }
    for (int i = 0; i < k; i++) {
        out << taus[i] << " " << iterations2[i] <<  endl;
    }
    for (int i = 0; i < k; i++) {
        out << taus[i] << " "  << iterations3[i] <<  endl;
    }
    out.close();
}


int main()
{
    int SystemNumber;
    int SizeTriangleMatrix;
    vector<vector<double>> TriangleMatrix{};
    vector<int> size{};
    vector<vector<vector<double>>> Matrix{};

    DataRead(SystemNumber, size, Matrix, "System.txt");


    WriteSimpleIterationAnswer(SystemNumber, size, Matrix, "SimpleIterationAnswer.txt");
    WriteSeidelAnswer(SystemNumber, size, Matrix, "SeidelAnswer.txt");
    //WriteJacobyAnswer(SystemNumber, size, Matrix, "JacobyAnswer.txt");

    //WriteRelaxationAnswer(SystemNumber, size, Matrix, "RelaxationAnswer.txt");


    //GenerateThreeDiagonal(201, "triangleMatrix.txt");
    //DataReadTriangleMatrix(SizeTriangleMatrix, TriangleMatrix, "triangleMatrix.txt");
    //WriteRelaxationTriangleAnswer(SizeTriangleMatrix, TriangleMatrix, "RelaxationTriangle.txt");
    //OmegaVsIteration(SystemNumber, size, Matrix, "OmegaVsIteration.txt", "WoframOmegaVsIteration.txt");
    //TauVsIteration(Matrix[0], "tau_vs_iter.txt", "wolfram_tau_vs_iter.txt", 500);

}