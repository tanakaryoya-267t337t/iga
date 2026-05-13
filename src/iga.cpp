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

vector<double> set_insert_knot(vector<double> &knotvector,vector<int> &knotspan,int p)
{
	vector<double> insert_knot;
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
	return insert_knot;
}


vector<double> knot_insertion(int p, int n){
	vector<double> knot = set_open_knot(p,n);
	vector<int> knotspan = set_knotspan(knot);
	vector<double> insert_knot = set_insert_knot(knot,knotspan,p);
	int N = n;
	vector<double> c(n*n,0.0);

	for(int i = 0; i < insert_knot.size(); i++){
		vector<double> new_knot;
		int k;
		for(int j = 0; j < knot.size()-1; j++){
			new_knot.push_back(knot.at(j));
			if(knot.at(j) <= insert_knot.at(i) && insert_knot.at(i) < knot.at(j+1)){
				new_knot.push_back(insert_knot.at(i));
				k = j;
			}
		}
        new_knot.push_back(knot.back());
		vector<double> alpha;
		for(int j = 0; j < N + 1; j++){
			if(j <= k - p){
				alpha.push_back(1.0);
			}
			else if(j <= k){
				alpha.push_back((insert_knot.at(i) - knot.at(j)) / (knot.at(j + p) - knot.at(j)));
			}
			else{
				alpha.push_back(0.0);
			}
		}
		vector<double> cT(N*(N+1),0.0);
		for(int j = 0; j < N ; j++){
			cT.at(j * N + j) = alpha.at(j);
			cT.at((j+1) * N + j) = 1.0-alpha.at(j+1);
		}
		vector<double> C(n*(N+1),0.0);
		if(i == 0){
			C = cT;
		}
		else {
			for(int j = 0; j < N + 1; j++){
				for(int k = 0; k < n; k++){
					for(int l = 0; l < N; l++){
						C.at(j*n+k) += cT.at(j*N+l)*c.at(l*n+k); 
					}
				}
			}
		}
		c.resize(n*(N+1));
		c = C;
        knot = new_knot;
		N++;
	}
    c = matT(c,n,N);
    return c;
}

void periodic_circle(int p,int n, vector<double>& cp){
    int ne = n - p;
    int nen = p + 1;
    vector<double> C(nen*nen,0.0);
    C.at(0) = 1.0/6.0;
    C.at(1) = 0.0;
    C.at(2) = 0.0;
    C.at(3) = 0.0;
    C.at(4) = 2.0/3.0;
    C.at(5) = 2.0/3.0;
    C.at(6) = 1.0/3.0;
    C.at(7) = 1.0/6.0;
    C.at(8) = 1.0/6.0;
    C.at(9) = 1.0/3.0;
    C.at(10) = 2.0/3.0;
    C.at(11) = 2.0/3.0;
    C.at(12) = 0.0;
    C.at(13) = 0.0;
    C.at(14) = 0.0;
    C.at(15) = 1.0/6.0;

    vector<double> R(ne*nen,0.0);

    int nt = 10;
    double dt = (double)1.0/nt;
    for(int it = 0; it < nt + 1; it++){
        double t = it * dt;
        vector<double> B = bernstein_basis_function(p,t);
        for(int i = 0; i < ne; i++){
            for(int j = 0; j < nen; j++){
                for(int k = 0; k < nen; k++){
                    R.at(i * nen + j) += C.at(j*nen+k)*B.at(k);
                }
            }
        }
    }
    double W = 0;
    vector<double> w(ne*nen,0.0);
    for(int i = 0; i < w.size(); i++){
        w.at(i) = 1.0;
    }

    double N = 0.0;
    
    for(int i = 0; i < ne; i++){
        for(int j = 0; j < nen; j++){
            N += w.at(i * nen + j) * R.at(i * nen + j);
        }
    }

    vector<double> p_x;
    vector<double> p_y;
    for(int i = 0; i < ne; i++){
        for(int j = 0; j < nen; j++){
            double px =0.0;
            double py =0.0;
            px = cp.at(i + j) * R.at(i * nen + j) / N;
            py = cp.at(n + i + j) * R.at(i * nen + j) / N;
            p_x.push_back(px);
            p_y.push_back(py);
        }
    }

    double r_ave = 0.0;
    for(int i = 0; i < p_x.size(); i++){
        r_ave += sqrt(p_x.at(i)*p_x.at(i) + p_y.at(i)*p_y.at(i)) / p_x.size();
    }

    for(int i = 0; i < n; i++){
        cp.at(i) /= r_ave;
        cp.at(n + i) /= r_ave;
    }
}

vector<double> periodic_circle_nurbs(int p, int n, double u){
    int ne = n;
    int nen = p + 1;
    vector<double> C(nen*nen,0.0);
    C.at(0) = 1.0/6.0;
    C.at(1) = 0.0;
    C.at(2) = 0.0;
    C.at(3) = 0.0;
    C.at(4) = 2.0/3.0;
    C.at(5) = 2.0/3.0;
    C.at(6) = 1.0/3.0;
    C.at(7) = 1.0/6.0;
    C.at(8) = 1.0/6.0;
    C.at(9) = 1.0/3.0;
    C.at(10) = 2.0/3.0;
    C.at(11) = 2.0/3.0;
    C.at(12) = 0.0;
    C.at(13) = 0.0;
    C.at(14) = 0.0;
    C.at(15) = 1.0/6.0;

    vector<double> R(nen,0.0);

    vector<double> B = bernstein_basis_function(p,u);
    for(int i = 0; i < nen; i++){
        for(int j = 0; j < nen; j++){
            R.at(i) += C.at(i * nen + j) * B.at(j);
        }
    }
    return R;
}

vector<double> bezier_element(int p, int n,double u,int i){
    vector<double> R(p+1,0.0);
    vector<double> C = knot_insertion(p,n);
    vector<double> knot = set_open_knot(p,n);
    vector<int> knotspan = set_knotspan(knot);
    vector<double> insert_knot = set_insert_knot(knot, knotspan, p);
    int m = insert_knot.size();
    int N = n+m;

    for(int j = 0; j < p + 1; j++){
            vector<double> B = bernstein_basis_function(p,u);
            for(int k = 0; k < p+1; k++){
                R.at(j) += C.at(i*(N+p)+j*N+k) * B.at(k);
        }
    }
    return R;
}

vector<double> nurbs_iga(int np, int px, int nx, int py, int ny, int pz, int nz, vector<double> & cp){
    #if 0 
    // periodic_circle(px,nx,cp);
    int ne = nx;
	vector<double> C(2 * np * ne, 0.0);
	vector<double> w(nx, 0.0);
	for (int i = 0; i < nx; i++)
	{
		w.at(i) = 1.0;
	}
    for(int ie = 0; ie < ne; ie++){
    	for (int i = 0; i < np; i++)
    	{
            double u = (double)1.0/(np-1)*i;
            vector<double> N = periodic_circle_nurbs(px, nx, u);
    		double R = 0.0;
    		for (int j = 0; j < px+1; j++)
    		{
                int id = (ie + j + nx - 1) % nx;
    			R += w.at(id)* N.at(j);
    		}
    		double cx = 0.0;
    		double cy = 0.0;
    		for (int j = 0; j < px + 1; j++)
    		{
                int id = (ie + j + nx - 1) % nx;
    			cx += w.at(id) * N.at(j) * cp.at(id) / R;
    			cy += w.at(id) * N.at(j) * cp.at(nx + id) / R;
    		}
    		C.at(ie*np+i) = cx;
    		C.at(np*ne + ie*np + i) = cy;
    	}
    }

    #elif 0 
    int nex = nx;
    int ney = ny -py;
    int ne = nex * ney;
    int nen = px + 1;
	vector<double> C(2 * ne * np * np, 0.0);
	vector<double> w(nx*ny, 0.0);
	for (int i = 0; i < ny; i++)
	{
        for(int j = 0; j < nx; j++){
            int id = i * nx + j;
            w.at(id)= 1.0;
        }
	}

    for(int iey = 0; iey < ney; iey++){
        for(int iex = 0; iex < nex; iex++){
            int ie = iex * iey;
	        for (int i = 0; i < np; i++)
	        {
               double uy = (double)1.0/(np-1)*i;
	        	vector<double> Ny = bezier_element(py,ny,uy,iey);
	        	for (int j = 0; j < np; j++)
	        	{
                  double ux = (double)1.0/(np-1)*j;
	        		vector<double> Nx = periodic_circle_nurbs(px,nx,ux);
	        		double R = 0.0;
	        		for (int a = 0; a < nen; a++)
	        		{
	        			for (int b = 0; b < nen; b++)
	        			{
                            int ix = (iex + b - 1 + nx) % nx;
                            int iy = iey + a;
                            int gid = iy * nx + ix;
	        				double N = Ny.at(a) * Nx.at(b);
	        				R += w.at(gid) * N;
	        			}
	        		}
	        		double cx = 0.0;
	        		double cy = 0.0;
	        		for (int a = 0; a < nen; a++)
	        		{
	        			for (int b = 0; b < nen; b++)
	        			{
                            int ix = (iex + b - 1 + nx) % nx;
                            int iy = iey + a;
                            int gid = iy * nx + ix;
	        				double N = Ny.at(a) * Nx.at(b);
	        				cx += w.at(gid) * N * cp.at(gid) / R;
	        				cy += w.at(gid) * N * cp.at(gid + nx * ny) / R;
	        			}
	        		}
	        		C.at(ie * np * np + i * np + j) = cx;
	        		C.at(ne * np * np + ie * np * np + i * np + j) = cy;
	        	}
	        }
        }
    }
    #else
    int nex = nx;
    int ney = ny - py;
    int nez = nz - pz;
    int ne = nex * ney * nez;
    int nen = px + 1;
	vector<double> C(3 * ne * np * np * np, 0.0);

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

    for(int iez = 0; iez < nez; iez++){
        for(int iey = 0; iey < ney; iey++){
            for(int iex = 0; iex < nex; iex++){
                int ie = iez * ney * nex + iey * nex + iex;
	            for (int i = 0; i < np; i++)
	            {
                   double uz = (double)1.0/(np-1)*i;
	            	vector<double> Nz = bezier_element(pz,nz,uz,iez);
	            	for (int j = 0; j < np; j++)
	            	{
                       double uy = 1.0/(np-1)*j;
	            	vector<double> Ny = bezier_element(py,ny,uy,iey);
	            		for(int k = 0; k < np; k++){
                           double ux = 1.0/(np-1)*k;
	            			vector<double> Nx = periodic_circle_nurbs(px,nx,ux);
	            			double R = 0.0;
	            			for (int a = 0; a < nen; a++)
	            			{
	            				for (int b = 0; b < nen; b++)
	            				{
	            					for(int c = 0; c < nen; c++){
                                        int ix = (iex + c - 1 + nx) % nx;
                                        int iy = iey + b;
                                        int iz = iez + a;
                                        int gid = iz * ny * nx + iy * nx + ix;
	            						double N = Nz.at(a) * Ny.at(b) * Nx.at(c);
	            						R += w.at(gid) * N;
	            					}
	            				}
	            			}
	            			double cx = 0.0;
	            			double cy = 0.0;
	            			double cz = 0.0;
	            			for (int a = 0; a < nen; a++)
	            			{
	            				for (int b = 0; b < nen; b++)
	            				{
	            					for(int c = 0; c < nen; c++){
                                        int ix = (iex + c - 1 + nx) % nx;
                                        int iy = iey + b;
                                        int iz = iez + c;
                                        int gid = iz * ny * nx + iy * nx + ix;
	            						double N = Nz.at(a) * Ny.at(b) * Nx.at(c);
	            						cx += w.at(gid) * N * cp.at(gid) / R;
	            						cy += w.at(gid) * N * cp.at(nx * ny * nz + gid) / R;
	            						cz += w.at(gid) * N * cp.at(2 * nx * ny * nz + gid) / R;
	            					}
	            				}
	            			}
	            			C.at(ie * np * np * np + i * np * np + j * np + k) = cx;
	            			C.at(ne * np * np * np + ie * np * np * np + i * np * np + j * np + k) = cy;
	            			C.at(2 * ne * np * np * np + ie * np * np * np +  i * np * np + j * np + k) = cz;
	            		}
	            	}
	            }
            }
        }
    }
    #endif
    return C;
}