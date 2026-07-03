#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "mydef.h"

using namespace std;

int main()
{
	int p = 3;
	int nx = 18;
	int ny = 5;
	int nz = 8;

    #ifdef CURVE
        int ny = 1;
        int nz = 1;
    #elif defined(SURFACE)
        int nz = 1;
    #endif
    int nn = nx * ny * nz;

    #ifdef MESENTERY
        int nx2 = 8;
        int ny2 = ny1;
        int nz2 = nz1;
        int nn2 = nx2 * ny2 * nz2;
        nx += nx2;
        ny += ny2;
        nz += nz2;
        nn += nn2;
    #endif




#ifdef CURVE
	vector<double> cp(2*nx,0.0);
	for(int i = 0; i < nx; i++){
		cp.at(i) = 1.0*cos(2*M_PI/(nx)*i);
		cp.at(nx+i) = 1.0*sin(2*M_PI/(nx)*i);
	}
#elif defined(SURFACE)
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
	vector<double> cp(3 * (nn), 0.0);
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
				cp.at(nn + id) = r * sin(theta);
				cp.at(2*nn + id) = l;
			}
		}
	}
    #ifdef MESENTERY
    for(int i = 0; i < nz2; i++){
        for(int j = 0; j < ny2; j++){
            for(int k = 0; k < nx2; k++){
                double y = 1.00 + 1.50/(nx2-1)*k;
                double x = 0.05 - 0.10/(ny2-1)*j;
                double s = (double)k/(nx2-1);
                double zc = 2.0;
                double half_width = 2.0*(1.0 -0.75*s);
                double t = (double)i/(nz2-1);
                double l = zc + half_width * (2.0*t -1.0);
                // double l = 4.0/(nz2-1)*i*(1.0/(k+1))+(2.0-2.0/(k+1));
                int id = i*nx2*ny2 + j*nx2 + k;
                cp.at(nn1 + id) = x;
                cp.at((nn1+nn2) + nn1 + id) = y;
                cp.at(2*(nn1+nn2) + nn1 + id) = l;
            }
        }
    }
    #endif

#endif


int np = 10;
int ne_bezier = nx1 * (ny1 - p) * (nz1 - p) + nx2 * (ny2 - p) * (nz2 - p);
int nn = np * np * np * ne_bezier;
int ne = (np-1) * (np-1) * (np - 1) * ne_bezier;
IGA iga(p, nx1, ny1, nz1, nx2, ny2, nz2, cp, np);
vector<double> C = iga.nurbs_iga();

string filename = "tube_with_messentary.vtk";
Output(np, nn, ne, C, filename);
return 0;
}