#pragma once
#include <vector>

class IGA {
    private:
    int p;
    int nx;
    int ny;
    int nz;
    std::vector<double> & cp;
    int np;

    public:
    std::vector<double> bernstein_basis_function(double u);
    std::vector<double> set_open_knot(int n);
    std::vector<double> bspline_basis_function(int n, double u, std::vector<double> knot);
    std::vector<double> bspline_curve();
    std::vector<double> bspline_surface();
    std::vector<double> bspline_volume();
    std::vector<double> nurbs_curve();
    std::vector<double> nurbs_surface();
    std::vector<double> nurbs_volume();
    std::vector<int> set_knotspan(std::vector<double> & knot);
    std::vector<double> set_insert_knot(std::vector<double> & knot, std::vector<int> & knotspan);
    std::vector<double> knot_insertion(int n);
    void periodic_circle(int n);
    std::vector<double> periodic_circle_nurbs(int n, double u);
    std::vector<double> bezier_element(int n, double u, int ie);
    std::vector<double> nurbs_iga();
    IGA (int p_, int nx_, int ny_, int nz_, std::vector<double> & cp_, int np_);
};

// std::vector<double> bernstein_basis_function(int p, double u);
// std::vector<double> set_open_knot(int p, int n);
// std::vector<double> bspline_basis_function(int p, int n, double u, std::vector<double> &knot);
// std::vector<double> bspline_curve(int p, int n, std::vector<double> & cp);
// std::vector<double> bspline_surface(int np,int p,int nx, int ny,std::vector<double> & cp);
// std::vector<double> bspline_volume(int np, int p,int nx, int ny, int nz, std::vector<double> &cp);
// std::vector<double> nurbs_curve(int np, int p, int n, std::vector<double> & cp);
// std::vector<double> nurbs_surface(int np, int p, int nx, int ny, std::vector<double> & cp);
// std::vector<double> nurbs_volume(int np, int p, int nx, int ny, int nz, std::vector<double> &cp);
// std::vector<int> set_knotspan(std::vector<double> &knotvector);
// std::vector<double> set_insert_knot(std::vector<double> &knotvector,std::vector<int> &knotspan, int p);
// std::vector<double> knot_insertion(int p, int n);
// void periodic_circle(int p, int n, std::vector<double>& cp);
// std::vector<double> periodic_circle_nurbs(int p, int n, double u);
// std::vector<double> bezier_element(int p, int n, double u);
// std::vector<double> nurbs_iga(int np, int p, int nx, int ny, int nz, std::vector<double> & cp);