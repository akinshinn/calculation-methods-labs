#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;



void triangularize(vector<vector<double>>& matrix);
vector<double> readf_vector(string file);
void print_matrix(const vector<vector<double>>& matrix);
vector<vector<double>> read_matrix();
void write_file(const vector<vector<double>>& matrix, string file);
vector<vector<double>> readf_matrix(string file_matrix, string file_vector);
int get_main_element_row(const vector<vector<double>>& matrix, int var_row);
void permutate_rows(vector<vector<double>>& matrix);
vector<double> reverse_course(const vector<vector<double>>& matrix);
vector<double> Gauss_method(vector<vector<double>>& matrix);
vector<double> check_ans(const vector<vector<double>>& matrix, const vector<double>& x);