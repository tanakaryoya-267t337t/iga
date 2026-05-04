#pragma once
#include <vector>
#include <string>

void Output(int np,int nn, int ne, std::vector<double>& x, std::string& filename);
void output_gnuplot(int n, int nn, std::vector<std::vector<double>> &B, std::string &gnuplotname);