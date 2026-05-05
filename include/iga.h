#pragma once
#include <vector>

std::vector<double> bernstein_basis_function(int p, double u);
std::vector<double> set_open_knot(int p, int a);
std::vector<double> bspline_basis_function(int p, int a, double u, std::vector<double> &k);
std::vector<double> bspline_curve(int p, int a, std::vector<double> & cp);
std::vector<double> bspline_surface(int np,int px,int py,int nx, int ny,std::vector<double> & cp);
std::vector<double> bspline_volume(int np, int px, int py, int pz, int nx, int ny, int nz, std::vector<double> &cp);
std::vector<double> nurbs_curve(int np, int px, int nx, std::vector<double> & cp);
std::vector<double> nurbs_surface(int np, int px, int nx, int py, int ny, std::vector<double> & cp);
std::vector<double> nurbs_volume(int np, int px, int nx, int py, int ny,int pz, int nz, std::vector<double> &cp);
std::vector<int> set_knotspan(std::vector<double> &knotvector);
std::vector<double> set_insert_knot(std::vector<double> &knotvector,std::vector<int> &knotspan, int p);
void knot_insert(int p, int a, std::vector<double> &knot, std::vector<double> &cp, double u, std::vector<double>& c,int count);
void knot_insertion(int p, int n);