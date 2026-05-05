#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include "output.h"
#include "iga.h"
#include "bicg.h"

using namespace std;

int main()
{
	int p = 3;
	int nx = 50;
	int ny = 5;
	int nz = 10;
	int first_n = 7;

#if 1 
	vector<double> cp(2*nx,0.0);
	for(int i = 0; i < nx; i++){
		cp.at(i) = 1.0*cos(2*M_PI/(nx-1)*i);
		cp.at(nx+i) = 1.0*sin(2*M_PI/(nx-1)*i);
	}
#elif 0 
	vector<double> cp(2 * nx * ny, 0.0);
	for (int i = 0; i < ny; i++)
	{
		cp.at(i * nx + 0) = (i * 0.025 + 0.40) * 1.0;
		cp.at(nx * ny + i * nx + 0) = (i * 0.025 + 0.40) * 0.0;
		cp.at(i * nx + 1) = (i * 0.025 + 0.40) * 1.0;
		cp.at(nx * ny + i * nx + 1) = (i * 0.025 + 0.40) * 1.0;
		cp.at(i * nx + 2) = (i * 0.025 + 0.40) * 0.0;
		cp.at(nx * ny + i * nx + 2) = (i * 0.025 + 0.40) * 1.0;
		cp.at(i * nx + 3) = (i * 0.025 + 0.40) * (-1.0);
		cp.at(nx * ny + i * nx + 3) = (i * 0.025 + 0.40) * 1.0;
		cp.at(i * nx + 4) = (i * 0.025 + 0.40) * (-1.0);
		cp.at(nx * ny + i * nx + 4) = (i * 0.025 + 0.40) * 0.0;
		cp.at(i * nx + 5) = (i * 0.025 + 0.40) * (-1.0);
		cp.at(nx * ny + i * nx + 5) = (i * 0.025 + 0.40) * (-1.0);
		cp.at(i * nx + 6) = (i * 0.025 + 0.40) * 0.0;
		cp.at(nx * ny + i * nx + 6) = (i * 0.025 + 0.40) * (-1.0);
		cp.at(i * nx + 7) = (i * 0.025 + 0.40) * 1.0;
		cp.at(nx * ny + i * nx + 7) = (i * 0.025 + 0.40) * (-1.0);
		cp.at(i * nx + 8) = (i * 0.025 + 0.40) * 1.0;
		cp.at(nx * ny + i * nx + 8) = (i * 0.025 + 0.40) * 0.0;
	}
#else
	vector<double> cp(3 * nx * ny * nz, 0.0);
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nx; k++)
			{
				int id = i * nx * ny + j * nx + k;
				cp.at(id) = (j * 0.10/(ny-1) + 0.40) * cos(2.0 * M_PI /(nx-1) * k);
				cp.at(nx * ny * nz + id) = (j * 0.10/(ny-1) + 0.40) * sin(2.0 * M_PI / (nx-1) * k);
				cp.at(2 * nx * ny * nz + id) = 4.0 / 9.0 * i;
			}
		}
	}
#endif

// vector<double> knot_x = set_open_knot(p, nx);
// vector<double> knot_y = set_open_knot(p, ny);

int np = 101;
int nn = np * np * np;
int ne = (np - 1) * (np - 1) * (np - 1);
vector<double> C = nurbs_volume(np, p, nx, p, ny, p, nz, cp);
vector<double> B;
vector<double> C(2 * nn, 0.0);
vector<double> cx;
vector<double> cy;
// 
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
// 
vector<vector<double>> b(n, vector<double>(nn, 0.0)); // plot
for (int i = 0; i < nn; i++)
{
double xi = 4.0 * i / (double)ne; // global parameter

int e = 0;
if (xi >= 0.0 && xi < 1.0)
e = 0;
else if (xi >= 1.0 && xi < 2.0)
e = 1;
else if (xi >= 2.0 && xi < 3.0)
e = 2;
else
e = 3;

double xi_left = e;
double xi_right = e + 1.0;
double u_local = (xi - xi_left) / (xi_right - xi_left);

B = bernstein(p, u_local);

double x = 0.0;
double y = 0.0;
vector<double> N(p + 1, 0.0);
for (int j = 0; j <= p; j++)
{
int global_id = e * p + j;
x += B.at(j) * cp.at(global_id);
y += B.at(j) * cp.at(n + global_id);
b.at(global_id).at(i) = B.at(j);
}
cx.push_back(x);
cy.push_back(y);
}
for (int i = 0; i < nn; i++)
{
C.at(i) = cx.at(i);
C.at(nn + i) = cy.at(i);
}

string filename = "nurbs_volume.vtk";
// string gnuplotname = "bspline_surface_basis_function.svg";
Output(np, nn, ne, C, filename);
// output_gnuplot(n, nn, b, gnuplotname);
return 0;
}