#pragma once
#include <vector>

double dot(std::vector<double>& a, std::vector<double>& b);
double product(std::vector<double>& A, std::vector<double>& b, int i);
double norm(std::vector<double>& a);
std::vector<double> matT(std::vector<double>& A, int col, int row);
std::vector<double> bicg(std::vector<double>& A, std::vector<double>& b, int N, bool bc);