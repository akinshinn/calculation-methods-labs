#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;
vector<float> readf_vector_float(string file);

vector<double> readf_vector_double(string file);



vector<vector<double>> readf_matrix_double(string file_matrix, string file_vector);

vector<vector<float>> readf_matrix_float(string file_matrix, string file_vector);

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
    T max = matrix[var_row][var_row];
    int max_row = var_row;
    int n = matrix.size();
    for (int i = var_row + 1; i < n; ++i) {
        if (matrix[i][var_row] > max) {
            max = matrix[i][var_row];
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
    permutate_rows(matrix);
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





template<typename T>
T Discrepancy(const vector<T>& vec1, const vector<T>& vec2) {
    T ans = 0;
    T suma = 0;
    if (vec1.size() == vec2.size()) {
        for (int i = 0; i < vec1.size(); i++) {
            suma = suma + (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
        }
        ans = sqrt(suma);
    }
    else {
        cout << "–азные размерности у векторов, дл€ которых ищем нев€зку";
    }
    return ans;
}


template<typename T>
T NormOne(const vector<vector<T>>& mas) {
    // последний столбец не учитываем тк у нас матрица A|B
    T suma = 0;
    T MaxSuma = -numeric_limits<T>::max();
    for (int j = 0; j < mas[0].size() - 1; j++) {
        suma = 0;
        for (int i = 0; i < mas.size(); i++) {
            suma += mas[i][j];
        }
        MaxSuma = max(MaxSuma, suma);

    }
    return MaxSuma;
}


template<typename T>
T NormInf(const vector<vector<T>>& mas) {
    // последний столбец не учитываем тк у нас матрица A|B
    T suma = 0;
    T MaxSuma = -numeric_limits<T>::max();
    for (int i = 0; i < mas.size(); i++) {
        suma = 0;
        for (int j = 0; j < mas[0].size() - 1; j++) {
            suma += mas[i][j];
        }
        MaxSuma = max(MaxSuma, suma);
    }
    return MaxSuma;
}


template<typename T>
T conditionality(const vector<vector<T>>& matrix, const  vector<vector<T>>& inv, T(*norm)(const vector<vector<T>>&) = NormOne) {
    return norm(matrix) * norm(inv);
}


template <typename T>
vector<vector<T>> transpose(const vector<vector<T>>& matrix) {
    int n = matrix.size();
    vector<vector<T>> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = vector<T>(n);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            res[i][j] = matrix[j][i];
            res[j][i] = matrix[i][j];
        }
    }
    return res;
}




// —читаем matrix расширенной 
template<typename T> 
vector<vector<T>> Gauss_inverse(vector<vector<T>>& matrix) {
    int n = matrix.size();
    vector<vector<T>> inv(n);
    vector<T> eye_col(n);


    for (int i = 0; i < n; ++i) {
        eye_col[i] = 1;
        if ((i - 1) >= 0) {
            eye_col[i - 1] = 0;
        }
        vector<vector<T>> temp_matrix = matrix_b_concat(matrix, eye_col);
        inv[i] = Gauss_method(temp_matrix);
        for (int k = 0; k < n; ++k) cout << inv[i][k] << " ";
        cout << endl;

    }
    inv = transpose(inv);


    return inv;
}



// ¬се матрицы считаютс€ расширенными, без вектора правой части все матрицы квадратные
template <typename T>
vector<vector<T>> matrix_prod(const vector<vector<T>>& m1, const vector<vector<T>>& m2) {
    int n = m1.size();
    vector<vector<T>> prod(n);
    for (int i = 0; i < n; ++i) {
        prod[i] = vector<T>(n);
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
            prod[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }


    return prod;
}


template <typename T>
void write_file(const vector<vector<T>>& matrix, const vector<T>& x, const vector<vector<T>>& inv, string file) {
    ofstream out;
    out.open(file);
    if (out.is_open()) {
        int n = matrix.size();
        out << n << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n + 1; ++j) {
                out << matrix[i][j] << " ";
            }
            out << endl;
        }
        out << endl;

        out << "solution" << endl;
        for (int i = 0; i < n; ++i) {
            out << x[i] << " ";
        }
        out << endl;

        out << "discrepancy" << endl;
        vector<T> b(n);
        for (int i = 0; i < n; ++i) b[i] = matrix[i][n];
        out << Discrepancy(x, b) << endl;

        out << "conditionality for norm one" << endl;
        out << conditionality(matrix, inv, NormOne)<<endl;

        out << "conditionality for norm inf" << endl;
        out << conditionality(matrix, inv, NormInf)<<endl;

        out.close();
    }
}