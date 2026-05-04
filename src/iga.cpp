#include <iostream>
#include <vector>
#include <cmath>
#include "output.h"
#include "iga.h"
#include "bicg.h"

using namespace std;

vector<double> bernstein_basis_function(int p, double u)
{
	vector<double> B(p + 1, 0.0);
	vector<double> Bn(p + 1, 0.0);
	B.at(0) = 1.0;
	for (int i = 0; i < p; i++)
	{
		Bn.at(0) = (1 - u) * B.at(0);
		for (int j = 1; j < p + 1; j++)
		{
			if (j > i + 1)
			{
				Bn.at(j) = 0.0;
			}
			else
			{
				Bn.at(j) = (1 - u) * B.at(j) + u * B.at(j - 1);
			}
		}
		for (int j = 0; j < p + 1; j++)
		{
			B.at(j) = Bn.at(j);
		}
	}
	return B;
}

vector<double> set_open_knot(int p, int a)
{
	vector<double> knot;
	double int_knot = 0.0;
	for (int i = 0; i <= p; i++)
	{
		knot.push_back(int_knot);
	}
	int num_int_knot = a + p + 1 - 2 * (p + 1);
	int_knot++;

	for (int i = 0; i < num_int_knot; i++)
	{
		knot.push_back(int_knot);
		int_knot++;
	}

	for (int i = 0; i <= p; i++)
	{
		knot.push_back(int_knot);
	}
	return knot;
}

vector<double> bspline_basis_function(int p, int a, double u, vector<double> &k)
{
	vector<double> N(a, 0.0);
	vector<double> Nn(a, 0.0);
	for (int i = 0; i < a; i++)
	{
		if (u >= k.at(i) && u < k.at(i + 1))
		{
			N.at(i) = 1.0;
		}
		else
		{
			N.at(i) = 0.0;
		}
	}
	if (u == k.back())
	{
		N.at(a - 1) = 1.0;
	}
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < a; j++)
		{
			double den1 = (k.at(j + i + 1) - k.at(j));
			double den2 = (k.at(j + i + 2) - k.at(j + 1));
			double turm1 = 0.0;
			double turm2 = 0.0;
			if (den1 != 0.0)
			{
				turm1 = (u - k.at(j)) / den1 * N.at(j);
			}
			if (j < a - 1 && den2 != 0.0)
			{
				turm2 = (k.at(j + i + 2) - u) / den2 * N.at(j + 1);
			}
			Nn.at(j) = turm1 + turm2;
		}
		for (int j = 0; j < a; j++)
		{
			N.at(j) = Nn.at(j);
		}
	}
	return N;
}

vector<double> bspline_curve(int p, int a, vector<double> &cp)
{
	int np = 101;
	vector<double> C(2 * np, 0.0);
	vector<double> knot = set_open_knot(p, a);
	double ks = knot.front();
	double ke = knot.back();
	for (int i = 0; i < np; i++)
	{
		double u = ks + (double)(ke - ks) * i / (np - 1);
		vector<double> N = bspline_basis_function(p, a, u, knot);
		double cx = 0.0;
		double cy = 0.0;
		for (int j = 0; j < a; j++)
		{
			cx += N.at(j) * cp.at(j);
			cy += N.at(j) * cp.at(a + j);
		}
		C.at(i) = cx;
		C.at(np + i) = cy;
	}
	return C;
}

vector<double> bspline_surface(int np, int px, int py, int nx, int ny, vector<double> &cp)
{
	vector<double> C(2 * np * np, 0.0);
	vector<double> knot_x = set_open_knot(px, nx);
	vector<double> knot_y = set_open_knot(py, ny);
	double ksx = knot_x.front();
	double kex = knot_x.back();
	double ksy = knot_y.front();
	double key = knot_y.back();
	for (int i = 0; i < np; i++)
	{
		double uy = ksy + (key - ksy) * i / (np - 1);
		vector<double> Ny = bspline_basis_function(py, ny, uy, knot_y);
		for (int j = 0; j < np; j++)
		{
			double ux = ksx + (kex - ksx) * j / (np - 1);
			vector<double> Nx = bspline_basis_function(px, nx, ux, knot_x);
			double cx = 0.0;
			double cy = 0.0;
			for (int a = 0; a < ny; a++)
			{
				for (int b = 0; b < nx; b++)
				{
					cx += Nx.at(b) * Ny.at(a) * cp.at(a * nx + b);
					cy += Nx.at(b) * Ny.at(a) * cp.at(nx * ny + a * nx + b);
				}
			}
			C.at(i * np + j) = cx;
			C.at(np * np + i * np + j) = cy;
		}
	}
	return C;
}

vector<double> bspline_volume(int np, int px, int py, int pz, int nx, int ny, int nz, vector<double> &cp)
{
	vector<double> C(3 * np * np * np, 0.0);
	vector<double> knot_x = set_open_knot(px, nx);
	vector<double> knot_y = set_open_knot(py, ny);
	vector<double> knot_z = set_open_knot(pz, nz);
	double ksx = knot_x.front();
	double kex = knot_x.back();
	double ksy = knot_y.front();
	double key = knot_y.back();
	double ksz = knot_z.front();
	double kez = knot_z.back();
	for (int i = 0; i < np; i++)
	{
		double uz = ksz + (kez - ksz) * i / (np - 1);
		vector<double> Nz = bspline_basis_function(pz, nz, uz, knot_z);
		for (int j = 0; j < np; j++)
		{
			double uy = ksy + (key - ksy) * j / (np - 1);
			vector<double> Ny = bspline_basis_function(py, ny, uy, knot_y);
			for (int k = 0; k < np; k++)
			{
				double ux = ksx + (kex - ksx) * k / (np - 1);
				vector<double> Nx = bspline_basis_function(px, nx, ux, knot_x);
				double cx = 0.0;
				double cy = 0.0;
				double cz = 0.0;
				for (int a = 0; a < nz; a++)
				{
					for (int b = 0; b < ny; b++)
					{
						for (int c = 0; c < nx; c++)
						{
							cx += Nx.at(c) * Ny.at(b) * Nz.at(a) * cp.at(a * nx * ny + b * nx + c);
							cy += Nx.at(c) * Ny.at(b) * Nz.at(a) * cp.at(nx * ny * nz + a * nx * ny + b * nx + c);
							cz += Nx.at(c) * Ny.at(b) * Nz.at(a) * cp.at(2 * nx * ny * nz + a * nx * ny + b * nx + c);
						}
					}
				}
				C.at(i * np * np + j * np + k) = cx;
				C.at(np * np * np + i * np * np + j * np + k) = cy;
				C.at(2 * np * np * np + i * np * np + j * np + k) = cz;
			}
		}
	}
	return C;
}

vector<double> nurbs_curve(int np, int px, int nx, vector<double> &cp)
{
	vector<double> C(2 * np, 0.0);
	vector<double> knot = set_open_knot(px, nx);
	vector<double> w(nx, 0.0);
	for (int i = 0; i < nx; i++)
	{
		w.at(i) = 1.0;
	}
	double ksx = knot.front();
	double kex = knot.back();
	for (int i = 0; i < np; i++)
	{
		double u = ksx + (double)(kex - ksx) * i / (np - 1);
		vector<double> N = bspline_basis_function(px, nx, u, knot);
		double R = 0.0;
		for (int j = 0; j < nx; j++)
		{
			R += w.at(j) * N.at(j);
		}
		double cx = 0.0;
		double cy = 0.0;
		for (int j = 0; j < nx; j++)
		{
			cx += w.at(j) * N.at(j) * cp.at(j) / R;
			cy += w.at(j) * N.at(j) * cp.at(nx + j) / R;
		}
		C.at(i) = cx;
		C.at(np + i) = cy;
	}
	return C;
}

vector<double> nurbs_surface(int np, int px, int nx, int py, int ny, vector<double> &cp)
{
	vector<double> C(2 * np * np, 0.0);
	vector<double> knot_x = set_open_knot(px, nx);
	vector<double> knot_y = set_open_knot(py, ny);

	vector<double> w(nx * ny, 0.0);
	for (int i = 0; i < ny; i++)
	{
		w.at(i * nx + 0) = 1.0;
		w.at(i * nx + 1) = 1.0 / sqrt(2.0);
		w.at(i * nx + 2) = 1.0;
		w.at(i * nx + 3) = 1.0 / sqrt(2.0);
		w.at(i * nx + 4) = 1.0;
		w.at(i * nx + 5) = 1.0 / sqrt(2.0);
		w.at(i * nx + 6) = 1.0;
		w.at(i * nx + 7) = 1.0 / sqrt(2.0);
		w.at(i * nx + 8) = 1.0;
	}

	double ksx = knot_x.front();
	double kex = knot_x.back();
	double ksy = knot_y.front();
	double key = knot_y.back();
	for (int i = 0; i < np; i++)
	{
		double uy = ksy + (double)(key - ksy) * i / (np - 1);
		vector<double> Ny = bspline_basis_function(py, ny, uy, knot_y);
		for (int j = 0; j < np; j++)
		{
			double ux = ksx + (double)(kex - ksx) * j / (np - 1);
			vector<double> Nx = bspline_basis_function(px, nx, ux, knot_x);
			double R = 0.0;
			vector<double> N(nx * ny, 0.0);
			for (int a = 0; a < ny; a++)
			{
				for (int b = 0; b < nx; b++)
				{
					int id = a * nx + b;
					N.at(id) = Ny.at(a) * Nx.at(b);
					R += w.at(id) * N.at(id);
				}
			}
			double cx = 0.0;
			double cy = 0.0;
			for (int a = 0; a < ny; a++)
			{
				for (int b = 0; b < nx; b++)
				{
					int id = a * nx + b;
					cx += w.at(id) * N.at(id) * cp.at(id) / R;
					cy += w.at(id) * N.at(id) * cp.at(nx * ny + id) / R;
				}
			}
			C.at(i * np + j) = cx;
			C.at(np * np + i * np + j) = cy;
		}
	}
	return C;
}

vector<double> nurbs_volume(int np, int px, int nx, int py, int ny,int pz, int nz, vector<double> &cp)
{
	vector<double> C(3 * np * np * np, 0.0);
	vector<double> knot_x = set_open_knot(px, nx);
	vector<double> knot_y = set_open_knot(py, ny);
	vector<double> knot_z = set_open_knot(pz, nz);

	vector<double> w(nx * ny * nz, 0.0);
	for (int i = 0; i < nz; i++)
	{
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nx; k++){
				int id = i * nx * ny + j * nx + k;
				w.at(id) = 1.0;
			}
		}
	}

	double ksx = knot_x.front();
	double kex = knot_x.back();
	double ksy = knot_y.front();
	double key = knot_y.back();
	double ksz = knot_z.front();
	double kez = knot_z.back();
	
	for (int i = 0; i < np; i++)
	{
		double uz = ksz + (double)(kez - ksz) * i / (np - 1);
		vector<double> Nz = bspline_basis_function(pz, nz, uz, knot_z);
		for (int j = 0; j < np; j++)
		{
		double uy = ksy + (double)(key - ksy) * j / (np - 1);
		vector<double> Ny = bspline_basis_function(py, ny, uy, knot_y);
			for(int k = 0; k < np; k++){
				double ux = ksx + (double)(kex - ksx) * k / (np - 1);
				vector<double> Nx = bspline_basis_function(px, nx, ux, knot_x);
				double R = 0.0;
				vector<double> N(nx * ny * nz, 0.0);
				for (int a = 0; a < nz; a++)
				{
					for (int b = 0; b < ny; b++)
					{
						for(int c = 0; c < nx; c++){
							int id = a * nx * ny + b * nx + c;
							N.at(id) = Nz.at(a) * Ny.at(b) * Nx.at(c);
							R += w.at(id) * N.at(id);
						}
					}
				}
				double cx = 0.0;
				double cy = 0.0;
				double cz = 0.0;
				for (int a = 0; a < nz; a++)
				{
					for (int b = 0; b < ny; b++)
					{
						for(int c = 0; c < nx; c++){
							int id = a * nx * ny + b *nx + c;
							cx += w.at(id) * N.at(id) * cp.at(id) / R;
							cy += w.at(id) * N.at(id) * cp.at(nx * ny * nz + id) / R;
							cz += w.at(id) * N.at(id) * cp.at(2 * nx * ny * nz + id) / R;
						}
					}
				}
				C.at(i * np * np + j * np + k) = cx;
				C.at(np * np * np + i * np * np + j * np + k) = cy;
				C.at(2 * np * np * np + i * np * np + j * np + k) = cz;
			}
		}
	}
	return C;
}

vector<int> set_knotspan(vector<double> &knotvector)
{
	vector<int> knotspan;
	double tmp = knotvector.front();
	for (int i = 0; i < knotvector.size(); i++)
	{
		if (tmp != knotvector.at(i))
		{
			knotspan.push_back(i);
			tmp = knotvector.at(i);
		}
	}
	return knotspan;
}

void set_insert_knot(
		vector<double> &knotvector,
		vector<int> &knotspan,
		vector<double> &insert_knot,
		int p)
{
	for (int i = 0; i < knotspan.size(); i++)
	{
		int a = knotspan.at(i);
		double s = knotvector.at(a);
		for (int j = 0; j < p; j++)
		{
			if (knotvector.at(a + j) != s)
			{
				insert_knot.push_back(s);
			}
		}
	}
}

void knot_insert(int p, int a, vector<double> &knot, vector<double> &cp, double u, vector<double> &c, int count)
{
	int first_n = 7;
	int n = knot.size();
	int ip = 0;
	vector<double> new_knot(n + 1, 0.0);
	for (int k = 0; k < n - 1; k++)
	{
		if (u >= knot.at(k) && u < knot.at(k + 1))
		{
			ip = k;
			new_knot.at(ip + 1) = u;
			int N = n - (ip + 1);
			for (int j = 0; j < n + 1; j++)
			{
				if (j <= ip)
				{
					new_knot.at(j) = knot.at(j);
				}
				else if (j > ip + 1)
				{
					new_knot.at(j) = knot.at(j - 1);
				}
			}
			break;
		}
	}

	vector<double> CP(2 * (a + 1), 0.0);
	vector<double> C((a + 1) * a, 0.0);

	double alpha = 0.0;

	for (int i = 0; i < a + 1; i++)
	{
		if (i <= ip - p)
		{
			CP.at(i) = cp.at(i);
			CP.at(a + 1 + i) = cp.at(a + i);
			alpha = 1.0;
		}
		else if (i > ip - p && i <= ip)
		{
			alpha = (double)(u - knot.at(i)) / (knot.at(i + p) - knot.at(i));
			CP.at(i) = alpha * cp.at(i) + (1.0 - alpha) * cp.at(i - 1);
			CP.at(a + 1 + i) = alpha * cp.at(a + i) + (1.0 - alpha) * cp.at(a + i - 1);
		}
		else
		{
			CP.at(i) = cp.at(i - 1);
			CP.at(a + 1 + i) = cp.at(a + i - 1);
			alpha = 0.0;
		}
		if (i < a)
		{
			C.at(i * (a + 1) + i) = alpha;
		}

		if (i > 0)
		{
			C.at((i - 1) * (a + 1) + i) = 1.0 - alpha;
		}
	}

	vector<double> CT = matT(C, a + 1, a);

	vector<double> c_new((a + 1) * (a - 1), 0.0);
	if (count == 1)
	{
		for (int i = 0; i < a + 1; i++)
		{
			for (int j = 0; j < first_n; j++)
			{
				for (int k = 0; k < a; k++)
				{
					c_new.at(i * first_n + j) += CT.at(i * a + k) * c.at(k * first_n + j);
				}
			}
		}
	}

	knot = new_knot;
	cp = CP;
	if (count == 0)
	{
		c = CT;
	}
	else if (count == 1)
	{
		c = c_new;
	}
}