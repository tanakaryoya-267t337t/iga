#include <vector>
#include <string>
#include <iostream>
#include "output.h"
#include "iga.h"
#include "bicg.h"

using namespace std;

int main()
{
	int p = 3;
	int n = 7;
	int first_n = 7;
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

	vector<int> knotspan = set_knotspan(knot);
	vector<double> insert_knot;
	set_insert_knot(knot, knotspan, insert_knot, p);

	int insert_knot_size = insert_knot.size();

	vector<double> c;
	int count = 0;
	for (int i = 0; i < insert_knot_size; i++)
	{
		knot_insert(p, n, knot, cp, insert_knot.at(i), c, count);
		n++;
		count = 1;
	}

	vector<double> cex = matT(c, first_n, n);

	vector<vector<double>> b(n, vector<double>(nn, 0.0)); // plot
	for (int i = 0; i < nn; i++)
	{
		double xi = 4.0 * i / (double)ne;  // global parameter

		int e = 0;
		if (xi >= 0.0 && xi < 1.0) e = 0;
		else if (xi >= 1.0 && xi < 2.0) e = 1;
		else if (xi >= 2.0 && xi < 3.0) e = 2;
		else e = 3;

		double xi_left = e;
		double xi_right = e + 1.0;
		double u_local = (xi - xi_left) / (xi_right - xi_left);

		B = bernstein(p, u_local);

		double x = 0.0;
		double y = 0.0;
		vector<double> N(p+1, 0.0);
		// for (int j = 0; j <= p; j++)
		// {
		// 	for (int k = 0; k <= p; k++)
		// 	{
		// 		N.at(j) += cex.at((e + j) * first_n + (e + k)) * B.at(k);
		// 	}
		// }
		for (int j = 0; j <= p ; j++)
		{
			int global_id = e*p + j;
			x += B.at(j) * cp.at(global_id);
			y += B.at(j) * cp.at(n + global_id);
			b.at(global_id).at(i) = N.at(j);
		}
		cx.push_back(x);
		cy.push_back(y);
	}
	for (int i = 0; i < nn; i++)
	{
		C.at(i) = cx.at(i);
		C.at(nn + i) = cy.at(i);
	}

	string filename = "bezeir_extraction_operator.vtk";
	string gnuplotname = "bezeir_extraction_operator_basis_function.svg";
	Output(nn, ne, C, filename);
	output_gnuplot(n, nn, b, gnuplotname);
	return 0;
}