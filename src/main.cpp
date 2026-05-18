#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "output.h"
#include "iga.h"
#include "bicg.h"

using namespace std;

int main()
{
	int p = 3;
	int nx1 = 10;
	int ny1 = 5;
	int nz1 = 8;
    int nn1 = nx1 * ny1 * nz1;

    int nx2 = 5;
    int ny2 = 8;
    int nz2 = nz1;
    int nn2 = nx2 * ny2 * nz2;

    int nx = nx1 + nx2;
    int ny = ny1 + ny2;
    int nz = nz1 + nz2;



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
	vector<double> cp(3 * (nn1+nn2), 0.0);
    double pi = acos(-1.0);
	for (int i = 0; i < nz1; i++)
	{
		for (int j = 0; j < ny1; j++)
		{
			for (int k = 0; k < nx1; k++)
			{
                double theta = 2.0*pi*k/nx1;
                double r = 0.80 + (1.0 - 0.80)/(ny1-1) * j;
                double l = 4.0/(nz1-1)*i;
				int id = i * nx1 * ny1 + j * nx1 + k;
				cp.at(id) = r * cos(theta);
				cp.at((nn1+nn2) + id) = r * sin(theta);
				cp.at(2*(nn1+nn2) + id) = l;
			}
		}
	}
    for(int i = 0; i < nz2; i++){
        for(int j = 0; j < ny2; j++){
            for(int k = 0; k < nx2; k++){
                double x = 0.10 - 0.20/(nx2-1)*k;
                double y = 1.0 + 1.50/(ny2-1)*j;
                double l = 4.0/(nz2-1)*i;
                int id = i*nx2*ny2 + j*nx2 + k;
                cp.at(nn1 + id) = x;
                cp.at((nn1+nn2) + nn1 + id) = y;
                cp.at(2*(nn1+nn2) + nn1 + id) = l;
            }
        }
    }

#endif


int np = 10;
int ne_bezier = nx1 * (ny1 - p) * (nz1 - p) + (nx2 - p) * (ny2 - p) * (nz2 - p);
int nn = np * np * np * ne_bezier;
int ne = (np-1) * (np-1) * (np - 1) * ne_bezier;
IGA iga(p, nx1, ny1, nz1, nx2, ny2, nz2, cp, np);
vector<double> C = iga.nurbs_iga();

string filename = "ien_set_messentary_tube.vtk";
Output(np, nn, ne, C, filename);
return 0;
}