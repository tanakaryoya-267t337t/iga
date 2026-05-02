#pragma once
#include <vector>

std::vector<double> bernstein(int p, double u);
std::vector<double> set_open_knot(int p, int a);
std::vector<double> bspline(int p, int a, double u, std::vector<double> &k);
std::vector<int> set_knotspan(std::vector<double> &knotvector);
void set_insert_knot(std::vector<double> &knotvector,std::vector<int> &knotspan,std::vector<double> &insert_knot,int p);
void knot_insert(int p, int a, std::vector<double> &knot, std::vector<double> &cp, double u, std::vector<double>& c,int count);