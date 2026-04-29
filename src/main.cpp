#include <vector>
#include <string>
#include <iostream>
#include "output.h"
#include "iga.h"

using namespace std;

int main()
{
	int p = 3;
	int n = 7;
	vector<double> cp(2 * n, 0.0);
	cp.at(0) = 1.0;
	cp.at(1) = 2.5;
	cp.at(2) = 3.4;
	cp.at(3) = 4.8;
	cp.at(4) = 5.2;
	cp.at(5) = 5.3;
	cp.at(6) = 6.0;
	cp.at(7) = 1.3;
	cp.at(8) = 2.6;
	cp.at(9) = 5.8;
	cp.at(10) = 2.6;
	cp.at(11) = 3.8;
	cp.at(12) = 6.9;
	cp.at(13) = 4.2;

	vector<double> knot = {0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4};

	int nn = 101;
	int ne = nn - 1;
	vector<double> B;
	vector<double> C(2 * nn, 0.0);
	vector<double> cx;
	vector<double> cy;


	vector<int>knotspan =  set_knotspan(knot);
	vector<double> insert_knot;
	set_insert_knot(knot,knotspan,insert_knot,p);

	int insert_knot_size = insert_knot.size();

	for(int i = 0; i < insert_knot_size; i++){
		knot_insert(p,n,knot,cp,insert_knot.at(i));
		n++;
	}


	vector<vector<double>> b(n,vector<double>(nn,0.0));  // plot
	for (int i = 0; i < nn; i++)
	{
		double u = 4.0 * i / (double)ne;
		B = bspline(p, n, u, knot);

		double x = 0.0;
		double y = 0.0;
		for (int a = 0; a < n; a++)
		{
			x += B.at(a) * cp.at(a);
			y += B.at(a) * cp.at(n + a);
			b.at(a).at(i) += B.at(a);
		}
		cx.push_back(x);
		cy.push_back(y);
	}
	for (int i = 0; i < nn; i++)
	{
		C.at(i) = cx.at(i);
		C.at(nn + i) = cy.at(i);
	}

	string filename = "knot_insertion.vtk";
	string gnuplotname = "knot_insertion_bspline_basis_function.svg";
	Output(nn, ne, C, filename);
	output_gnuplot(n, nn, b, gnuplotname);
	return 0;
}