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
	int nx = 10;
	int ny = 5;
	int nz = 8;
	int first_n = 7;

#if 0 
	vector<double> cp(2*nx,0.0);
	for(int i = 0; i < nx; i++){
		cp.at(i) = 1.0*cos(2*M_PI/(nx)*i);
		cp.at(nx+i) = 1.0*sin(2*M_PI/(nx)*i);
	}
#elif 0 
	vector<double> cp(2 * nx * ny, 0.0);
    double pi = acos(-1.0);
    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            double theta = 2.0*pi*j/nx;
            double r =0.80 + (1.0 - 0.80)/(ny-1) * i;
            int id = i * nx + j;
            cp.at(id) = r*cos(theta);
            cp.at(nx*ny+id) = r*sin(theta);
        }
    }
#else
	vector<double> cp(3 * nx * ny * nz, 0.0);
    double pi = acos(-1.0);
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nx; k++)
			{
                double theta = 2.0*pi*k/nx;
                double r = 0.80 + (1.0 - 0.80)/(ny-1) * j;
                double l = 4.0/(nz-1)*i;
				int id = i * nx * ny + j * nx + k;
				cp.at(id) = r * cos(theta);
				cp.at(nx * ny * nz + id) = r * sin(theta);
				cp.at(2 * nx * ny * nz + id) = l;
			}
		}
	}
#endif


int np = 21;
int ne_bezier = nx*(ny-p)*(nz-p);
int nn = np * np * np * ne_bezier;
int ne = (np-1) * (np-1) * (np - 1) * ne_bezier;
// int nn = np ;
// int ne = (np - 1);
vector<double> C = nurbs_iga(np, p, nx, p, ny, p, nz, cp);

string filename = "periodic_tube.vtk";
Output(np, nn, ne, C, filename);
// string gnuplotname = "bspline_surface_basis_function.svg";
// output_gnuplot(n, nn, b, gnuplotname);
return 0;
}