#pragma once
#include <vector>

class IGA {
    private:
    int p;
    int ncp1;
    int ncp2;
    int ncp3;
    int mcp1;
    int mcp2;
    int mcp3;
    int np;
    int nnt;
    int nnm;
    int nn;
    int net;
    int nem;
    int ne;
    int nen;
    std::vector<double> & cp;
    std::vector<int> ien;
    std::vector<int> ien_offset;

    public:
    std::vector<double> bernstein_basis_function(double u);
    std::vector<double> set_open_knot(int n);
    std::vector<double> bspline_basis_function(int n, double u, std::vector<double> knot);
    std::vector<double> bspline_curve(int ncp);
    std::vector<double> bspline_surface(int ncp1, int ncp2);
    std::vector<double> bspline_volume(int ncp1, int ncp2, int ncp3);
    std::vector<double> nurbs_curve(int ncp1);
    std::vector<double> nurbs_surface(int ncp1, int ncp2);
    std::vector<double> nurbs_volume(int ncp1, int ncp2, int ncp3);
    std::vector<int> set_knotspan(std::vector<double> & knot);
    std::vector<double> set_insert_knot(std::vector<double> & knot, std::vector<int> & knotspan);
    std::vector<std::vector<double>> knot_insertion(int n);
    void periodic_tube(int ncp);
    std::vector<double> periodic_circle_nurbs(int n, double u);
    std::vector<double> bezier_element(int n, double u, int ie, std::vector<std::vector<double>> & C);
    std::vector<double> connect_patch_basis_x(double tx, int iex, int ix,std::vector<std::vector<double>> & Cx);
    std::vector<double> nurbs_basis(int ie, bool tube);
    std::vector<double> nurbs_iga();
    IGA (int p_, int ncp1_, int ncp2_, int ncp3_, int mcp1_, int mcp2_, int mcp3_, std::vector<double> & cp_, int np_);
};
