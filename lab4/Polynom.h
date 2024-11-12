#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
const double epsilon = 0;

class Polynom {
	public: 
		int n;
		vector<double> coefs;
		void copy(const Polynom& P) {
			n = P.n;
			coefs = P.coefs;
		}

		Polynom() {
			n = 0;
			coefs.push_back(0);
		}
		Polynom(int new_n) {
			n = new_n;
			for (int i = 0; i < n + 1; i++) {
				coefs.push_back(0);
			}
		}
		Polynom(const vector<double> new_coefs) {
			coefs = new_coefs; 
			n = new_coefs.size() - 1;
		}
		Polynom(const Polynom& P) {
			this->copy(P);
		}
		void print() {
			for (int i = 0; i < n; i++) {
				if (coefs[i] != 1) {
					if (n - i == 1) {
						cout << coefs[i] << "*" << "x" << " + ";
					}
					else {
						cout << coefs[i] << "*" << "x^" << n - i << " + ";
					}
				}
				else {
					if (n - i == 1) {
						cout <<  "x" << " + ";
					}
					else {
						cout << "x^" << n - i << " + ";
					}
				}
					
			}
			cout << coefs[n] << endl;
		}
		void write_file(string file) {
			ofstream f;
			f.open(file);
			if (f.is_open()) {
				for (int i = 0; i < n + 1; ++i) {
					f << coefs[i] << " ";

				}
				f.close();
			}
			
		}
		void multiply(const Polynom& P2) {
			Polynom res(n + P2.n);
			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < P2.n + 1; j++) {
					res.coefs[i + j] += coefs[i] * P2.coefs[j];
				}
			}
			copy(res);
		}
		void add(const Polynom& P2) {
			int max_n = max(n, P2.n);
			int min_n = n + P2.n - max_n;
			Polynom& minP = *this, maxP = P2;
			Polynom res(maxP);

			if (max_n != min_n) {
				if (max_n == n){
					maxP = *this;
					minP = P2;
				}
				else {
					maxP = P2;
					minP = *this;
				}

				for (int i = min_n; i >= 0; --i) {
					res.coefs[i + (max_n - min_n)] += minP.coefs[i];
				}
			}
			else {
				for (int i = 0; i < max_n + 1; i++) {
					res.coefs[i] += minP.coefs[i];
				}
			}
			
			copy(res);
		}
		Polynom& operator+=(const Polynom& p2) {
			add(p2);
			return *this;
		}
		Polynom& operator*=(const Polynom& p2) {
			multiply(p2);
			return *this;
		}
		Polynom& operator*=(const double& a) {
			for (int i = 0; i < n + 1; i++) {
				coefs[i] *= a;
			}
			return *this;
		}
		Polynom& operator=(const Polynom& p2) {
			copy(p2);
			return *this;
		}

		double getValue(double x) {
			double res = 0;
			for (int i = 0; i < n + 1; i++) {
				res += coefs[i] * pow(x, n - i);
			}
			return res;
		}
};