#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;



void get_U(vector<vector<double>>* matrix);

void print_matrix(const vector<vector<double>>& matrix);
vector<vector<double>> read_matrix();
void write_file(const vector<vector<double>>& matrix, string file);
vector<vector<double>> readf_matrix(string file);