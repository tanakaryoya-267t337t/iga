#pragma once
#include <vector>

class IGA {
    private:
    int p;
    int nx;
    int ny;
    int nz;
    int np;
    int nn;
    int ne;
    int nen;
    std::vector<double> & cp;
    std::vector<int> ien;
    std::vector<int> ien_offset;

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
    std::vector<std::vector<double>> knot_insertion(int n);
    void periodic_circle(int n);
    std::vector<double> periodic_circle_nurbs(int n, double u);
    std::vector<double> bezier_element(int n, double u, int ie);
    std::vector<double> nurbs_basis(int ie);
    std::vector<double> nurbs_iga();
    IGA (int p_, int nx_, int ny_, int nz_, std::vector<double> & cp_, int np_);
};
