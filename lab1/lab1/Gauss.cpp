#include "Gauss.h"
//#include "QR.cpp"

vector<vector<double>> read_matrix() {
	int n;
	cin >> n;
	vector<vector<double>> matrix;

	for (int i = 0; i < n; ++i) {
		matrix.emplace_back(vector<double>(n));
		for (int j = 0; j < n; ++j) {
			cin >> matrix[i][j];
		}
	}
	return matrix;
}


vector<vector<float>> readf_matrix_float(string file_matrix) {
	ifstream in;
	in.open(file_matrix);
	int n;
	vector<vector<float>> matrix;
	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			matrix.emplace_back(vector<float>(n + 1));
			for (int j = 0; j < n + 1; ++j) {
				in >> matrix[i][j];
			}
		}
		in.close();
	}
	return matrix;
}









//vector<double> readf_vector_double(string file) {
//	ifstream in;
//	in.open(file);
//	int n;
//	vector<double> b;
//	double temp;
//	if (in.is_open()) {
//		in >> n;
//		for (int i = 0; i < n; ++i) {
//			in >> temp;
//			b.emplace_back(temp);
//		}
//		in.close();
//	}
//	return b;
//}




vector<vector<double>> readf_matrix_double(string file_matrix) {
	//vector<double> b = readf_vector_double(file_vector);
	ifstream in;
	in.open(file_matrix);
	int n;
	vector<vector<double>> matrix;

	if (in.is_open()) {
		in >> n;
		for (int i = 0; i < n; ++i) {
			matrix.emplace_back(vector<double>(n + 1));
			for (int j = 0; j < n+1; ++j) {
				in >> matrix[i][j];
			}
			//matrix[i][n] = b[i];
		}
		in.close();
	}
	return matrix;
}




pair<vector<float>, int> QRFloat(int size, vector<vector<float>>& Matrix) {
    // Прямой ход
    float c, s;
    float tmp1, tmp2;
    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        for (int j = i + 1; j < size; j++) {
            // j - index уравнения в котором будем устанять x_i
            if (Matrix[j][i] == 0) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);

            for (int k = 0; k < size + 1; k++) {
                // size + 1 потому что Matrix состоит из n+1 столбца Matrix = A | b
                tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
                tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
                Matrix[i][k] = tmp1;
                Matrix[j][k] = tmp2;
            }


        }
    }
    // Проверка
    float product = 1;
    for (int i = 0; i < size; i++) {
        product *= Matrix[i][i];
    }

    vector<float> ans(size);
    int flag = 1;
    if (product != 0) {
        // Обратный код
        int n = size - 1;
        float suma = 0;
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


pair<vector<float>, int> QRFloatDecomposition(int size, vector<vector<float>>& Matrix, vector<vector<float>>& TMatrix) {
    // Прямой ход
    float c, s;
    float tmp1, tmp2, tmp3, tmp4;
    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        for (int j = i + 1; j < size; j++) {
            // j - index уравнения в котором будем устанять x_i
            if (Matrix[j][i] == 0) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);

            for (int k = 0; k < size + 1; k++) {
                // size + 1 потому что Matrix состоит из n+1 столбца Matrix = A | b
                tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
                tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;

                Matrix[i][k] = tmp1;
                Matrix[j][k] = tmp2;
            }

            // пересчет iой строки

            for (int k = 0; k < size; k++) {
                tmp3 = c * TMatrix[i][k] + s * TMatrix[j][k];
                tmp4 = -s * TMatrix[i][k] + c * TMatrix[j][k];
                TMatrix[i][k] = tmp3;
                TMatrix[j][k] = tmp4;
            }
        }
    }
    // Проверка
    float product = 1;
    for (int i = 0; i < size; i++) {
        product *= Matrix[i][i];
    }

    vector<float> ans(size);
    int flag = 1;
    if (product != 0) {
        // Обратный код
        int n = size - 1;
        float suma = 0;
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

pair<vector<double>, int> QRDouble(int size, vector<vector<double>>& Matrix) {
    // Прямой ход
    double c, s;
    double tmp1, tmp2;
    for (int i = 0; i < size - 1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        for (int j = i + 1; j < size; j++) {
            // j - index уравнения в котором будем устанять x_i
            if (Matrix[j][i] == 0) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);

            for (int k = 0; k < size + 1; k++) {
                // size + 1 потому что Matrix состоит из n+1 столбца Matrix = A | b
                tmp1 = Matrix[i][k] * c + Matrix[j][k] * s;
                tmp2 = Matrix[j][k] * c - Matrix[i][k] * s;
                Matrix[i][k] = tmp1;
                Matrix[j][k] = tmp2;
            }
        }
    }
    // Проверка
    double product = 1;
    for (int i = 0; i < size; i++) {
        product *= Matrix[i][i];
    }

    vector<double> ans(size);
    int flag = 1;
    if (product != 0) {
        // Обратный код
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

pair<vector<double>, int> QRDoubleDecomposition(int size, vector<vector<double>>& Matrix, vector<vector<double>>& TMatrix) {
    // Прямой ход
    double c, s;
    double tmp1, tmp2, tmp3, tmp4;
    cout << "Decomposition" << endl;
    for (int i = 0; i < size-1; i++) {
        // i - index переменной которой будем устранять, если i = 2 устраняем x2 во всех уранениях
        for (int j = i + 1; j < size; j++) {
            // j - index уравнения в котором будем устанять x_i
            if (Matrix[j][i] < 1e-15) continue;
            c = Matrix[i][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);
            s = Matrix[j][i] / sqrt(Matrix[i][i] * Matrix[i][i] + Matrix[j][i] * Matrix[j][i]);

            for (int k = 0; k < size + 1; k++) {
                // size + 1 потому что Matrix состоит из n+1 столбца Matrix = A | b
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
            cout << "i = " << i << " j = " << j << endl;
            cout << "c = " << c << " s = " << s << endl;
            print_matrix(TMatrix);
        }
    }
    // Проверка
    double product = 1;
    for (int i = 0; i < size; i++) {
        product *= Matrix[i][i];
    }

    vector<double> ans(size);
    int flag = 1;
    if (product != 0) {
        // Обратный код
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

vector<vector<float>> InverseFloat(vector<vector<float>> Matrix) {
    int size = Matrix.size();
    vector<vector<float>> InverseMatrix{}, CopyMatrix{}, resMatrix{};

    for (int i = 0; i < size; i++) {
        Matrix[i][size] = 0;
    }
    for (int i = 0; i < size; i++) {
        Matrix[i][size] = 1;
        if (i != 0) Matrix[i - 1][size] = 0;
        CopyMatrix = Matrix;


        pair<vector<float>, int> ans = QRFloat(size, CopyMatrix);
        InverseMatrix.push_back(ans.first);
    }

    resMatrix.resize(size, vector<float>(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            resMatrix[i][j] = InverseMatrix[j][i];
        }
    }
    return resMatrix;
}


vector<vector<double>> InverseDouble(vector<vector<double>> Matrix) {
    int size = Matrix.size();
    vector<vector<double>> InverseMatrix{}, CopyMatrix{}, resMatrix{};

    for (int i = 0; i < size; i++) {
        Matrix[i][size] = 0;
    }
    for (int i = 0; i < size; i++) {
        Matrix[i][size] = 1;
        if (i != 0) Matrix[i - 1][size] = 0;
        CopyMatrix = Matrix;


        pair<vector<double>, int> ans = QRDouble(size, CopyMatrix);
        InverseMatrix.push_back(ans.first);

    }
    resMatrix.resize(size, vector<double>(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            resMatrix[i][j] = InverseMatrix[j][i];
        }
    }
    return resMatrix;
}




void FloatTest(string filename) {
    int SystemNumber = 0;
    vector<int>  size{};
    vector<vector<vector<float>>> Matrix{}, InverseMatrix{}, CopyMatrix{}, CheckMatrix{}, RMatrix{}, TMatrix{};
    DataRead(SystemNumber, size, Matrix, filename);
    CopyMatrix = Matrix;


    vector<pair<vector<float>, int>> OutPutData = {};
    vector<vector<float>> TIdentityMatrix{};
    pair<vector<float>, int> vec;
    vector<float> CalculatedRightParts{}, RealRightParts{}, TrueRightParts{}, OutragedRigthPart{}, OutragedSol{}, vecDiscrepancy{}, ConditionalityOne{}, ConditionalityInfty{}, vecAssessmentConditionalityOne{}, vecAssessmentConditionalityInfty{};
    for (int k = 0; k < SystemNumber; k++) {
        TIdentityMatrix.resize(size[k], vector<float>(size[k], 0));
        for (int u = 0; u < size[k]; u++) {
            TIdentityMatrix[u][u] = 1;
        }

        vec = QRFloatDecomposition(size[k], Matrix[k], TIdentityMatrix);

        OutPutData.push_back(vec);


        RealRightParts = {};
        if (vec.second == 1) {
            RMatrix.push_back(Matrix[k]);
            TMatrix.push_back(TIdentityMatrix);

            InverseMatrix.push_back(InverseFloat(CopyMatrix[k]));
            CheckMatrix.push_back(MultiplyMatrix(InverseMatrix[k], CopyMatrix[k]));

            //невязка
            CalculatedRightParts = SearchRightParts(vec.first, CopyMatrix[k]);
            for (int ind = 0; ind < Matrix[k].size(); ind++) {
                RealRightParts.push_back(CopyMatrix[k][ind][size[k]]);
            }
            vecDiscrepancy.push_back(Discrepancy(RealRightParts, CalculatedRightParts));

            //Обусловленность Можно выбрать одну из норм NormOne или NormIfty
            ConditionalityOne.push_back(SearchConditionality(CopyMatrix[k], InverseMatrix[k], NormOneMatrix));
            ConditionalityInfty.push_back(SearchConditionality(CopyMatrix[k], InverseMatrix[k], NormInfMatrix));

            //Возмущаем систему
            TrueRightParts = {};
            OutragedRigthPart = {};
            OutragedSol = {};
            for (int i = 0; i < size[k]; i++) {
                TrueRightParts.push_back(CopyMatrix[k][i][size[k]]);

                if (i % 2 == 0) {
                    CopyMatrix[k][i][size[k]] += outrageFloat;
                }

                OutragedRigthPart.push_back(CopyMatrix[k][i][size[k]]);
            }

            OutragedSol = QRFloat(size[k], CopyMatrix[k]).first;

            vecAssessmentConditionalityOne.push_back(AssessmentConditionality(vec.first, OutragedSol, TrueRightParts, OutragedRigthPart, NormOneVec));
            vecAssessmentConditionalityInfty.push_back(AssessmentConditionality(vec.first, OutragedSol, TrueRightParts, OutragedRigthPart, NormInfVec));
        }
        else {
            InverseMatrix.push_back({});
            vecDiscrepancy.push_back(flagOfNoSolution);
            ConditionalityOne.push_back(flagOfNoSolution);
            ConditionalityInfty.push_back(flagOfNoSolution);
            vecAssessmentConditionalityOne.push_back(flagOfNoSolution);
            vecAssessmentConditionalityInfty.push_back(flagOfNoSolution);
            RMatrix.push_back({});
            TMatrix.push_back({});
        }



    }

    WriteAnswer(OutPutData, "Answer.txt");

    WriteDiscrepancy(vecDiscrepancy, "Discrepancy.txt");

    WriteConditionality(ConditionalityOne, ConditionalityInfty, "Conditionality.txt");

    WriteCheckPoint(CheckMatrix, "CheckPoint.txt");

    WriteAssessmentConditionality(vecAssessmentConditionalityOne, vecAssessmentConditionalityInfty, "AssessmentConditionality.txt");

    WriteQRMatrix(RMatrix, TMatrix, "QRDecomposition.txt");
}



void DoubleTest(string filename) {
    int SystemNumber = 0;
    vector<int>  size{};
    vector<vector<vector<double>>> Matrix{}, InverseMatrix{}, CopyMatrix{}, CheckMatrix{}, RMatrix{}, TMatrix{};
    DataRead(SystemNumber, size, Matrix, filename);

    CopyMatrix = Matrix;

    vector<pair<vector<double>, int>> OutPutData = {};
    vector<vector<double>> TIdentityMatrix{};
    pair<vector<double>, int> vec;
    vector<double>  CalculatedRightParts{}, RealRightParts{}, TrueRightParts{}, OutragedRigthPart{}, OutragedSol{}, vecDiscrepancy{}, ConditionalityOne{}, ConditionalityInfty{}, vecAssessmentConditionalityOne{}, vecAssessmentConditionalityInfty{};

    for (int k = 0; k < SystemNumber; k++) {
        TIdentityMatrix.resize(size[k], vector<double>(size[k], 0));
        for (int u = 0; u < size[k]; u++) {
            TIdentityMatrix[u][u] = 1;
        }
        vec = QRDoubleDecomposition(size[k], Matrix[k], TIdentityMatrix);

        OutPutData.push_back(vec);

        RealRightParts = {};
        if (vec.second == 1) {
            RMatrix.push_back(Matrix[k]);
            TMatrix.push_back(TIdentityMatrix);


            InverseMatrix.push_back(InverseDouble(CopyMatrix[k]));
            CheckMatrix.push_back(MultiplyMatrix(InverseMatrix[k], CopyMatrix[k]));

            //невязка
            CalculatedRightParts = SearchRightParts(vec.first, CopyMatrix[k]);
            for (int ind = 0; ind < Matrix[k].size(); ind++) {
                RealRightParts.push_back(CopyMatrix[k][ind][size[k]]);
            }
            vecDiscrepancy.push_back(Discrepancy(RealRightParts, CalculatedRightParts));

            //Обусловленность Можно выбрать одну из норм NormOne или NormIfty
            ConditionalityOne.push_back(SearchConditionality(CopyMatrix[k], InverseMatrix[k], NormOneMatrix));
            ConditionalityInfty.push_back(SearchConditionality(CopyMatrix[k], InverseMatrix[k], NormInfMatrix));

            //Возмущаем систему
            TrueRightParts = {};
            OutragedRigthPart = {};
            OutragedSol = {};
            for (int i = 0; i < size[k]; i++) {
                TrueRightParts.push_back(CopyMatrix[k][i][size[k]]);

                if (i % 2 == 0) {
                    CopyMatrix[k][i][size[k]] += outrageFloat;
                }

                OutragedRigthPart.push_back(CopyMatrix[k][i][size[k]]);
            }

            OutragedSol = QRDouble(size[k], CopyMatrix[k]).first;

            vecAssessmentConditionalityOne.push_back(AssessmentConditionality(vec.first, OutragedSol, TrueRightParts, OutragedRigthPart, NormOneVec));
            vecAssessmentConditionalityInfty.push_back(AssessmentConditionality(vec.first, OutragedSol, TrueRightParts, OutragedRigthPart, NormInfVec));
        }
        else {
            InverseMatrix.push_back({});
            vecDiscrepancy.push_back(flagOfNoSolution);
            ConditionalityOne.push_back(flagOfNoSolution);
            ConditionalityInfty.push_back(flagOfNoSolution);
            vecAssessmentConditionalityOne.push_back(flagOfNoSolution);
            vecAssessmentConditionalityInfty.push_back(flagOfNoSolution);
            RMatrix.push_back({});
            TMatrix.push_back({});
        }
    }
    WriteAnswer(OutPutData, "Answer.txt");

    WriteDiscrepancy(vecDiscrepancy, "Discrepancy.txt");

    WriteConditionality(ConditionalityOne, ConditionalityInfty, "Conditionality.txt");

    WriteCheckPoint(CheckMatrix, "CheckPoint.txt");

    WriteAssessmentConditionality(vecAssessmentConditionalityOne, vecAssessmentConditionalityInfty, "AssessmentConditionality.txt");

    WriteQRMatrix(RMatrix, TMatrix, "QRDecomposition.txt");
}