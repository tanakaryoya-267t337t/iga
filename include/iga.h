#pragma once
#include <vector>

class IGA {
    private:
    int p;
    int nx1;
    int ny1;
    int nz1;
    int nx2;
    int ny2;
    int nz2;
    int np;
    int nn1;
    int nn2;
    int nn;
    int ne1;
    int ne2;
    int ne;
    int nen;
    std::vector<double> & cp;
    std::vector<int> ien;
    std::vector<int> ien_offset;

    public:
    std::vector<double> bernstein_basis_function(double u);
    std::vector<double> set_open_knot(int n);
    std::vector<double> bspline_basis_function(int n, double u, std::vector<double> knot);
    std::vector<double> bspline_curve(int nx);
    std::vector<double> bspline_surface(int nx, int ny);
    std::vector<double> bspline_volume(int nx, int ny, int nz);
    std::vector<double> nurbs_curve(int nx);
    std::vector<double> nurbs_surface(int nx, int ny);
    std::vector<double> nurbs_volume(int nx, int ny, int nz);
    std::vector<int> set_knotspan(std::vector<double> & knot);
    std::vector<double> set_insert_knot(std::vector<double> & knot, std::vector<int> & knotspan);
    std::vector<std::vector<double>> knot_insertion(int n);
    void periodic_tube(int nx);
    std::vector<double> periodic_circle_nurbs(int n, double u);
    std::vector<double> bezier_element(int n, double u, int ie, std::vector<std::vector<double>> & C);
    std::vector<double> connect_patch_basis_x(double tx, int iex, int ix,std::vector<std::vector<double>> & Cx);
    std::vector<double> nurbs_basis(int ie, bool tube);
    std::vector<double> nurbs_iga();
    IGA (int p_, int nx1_, int ny1_, int nz1_, int nx2_, int ny2_, int nz2_, std::vector<double> & cp_, int np_);
};
