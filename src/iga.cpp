#include <iostream>
#include <vector>
#include <cmath>
#include "output.h"
#include "iga.h"
#include "bicg.h"

using namespace std;

IGA::IGA (int p_, int nx1_, int ny1_, int nz1_,int nx2_, int ny2_, int nz2_, vector<double> & cp_, int np_):p(p_),nx1(nx1_),ny1(ny1_),nz1(nz1_),nx2(nx2_),ny2(ny2_),nz2(nz2_),cp(cp_),np(np_){
    this->nn = nx1_*ny1_*nz1_ + nx2_*ny2_*nz2_;
    this->nen = p + 1;

    int nex1 = nx1_;
    int ney1 = ny1_ - p;
    int nez1 = nz1_ - p;

    int nex2 = nx2_- p;
    int ney2 = ny2_ - p;
    int nez2 = nz2_ - p;

    this->nn1 = nx1_ * ny1_ * nz1_;
    this->nn2 = nx2_ * ny2_ * nz2_;


    this->ne1 = nex1 * ney1 * nez1;
    this->ne2 = nex2 * ney2 * nez2;
    this->ne = ne1 + ne2;


    for(int iez = 0; iez < nez1; iez++){
        for(int iey = 0; iey < ney1; iey++){
            for(int iex = 0; iex < nex1; iex++){
                for(int i = 0; i <nen; i++){
                    for(int j = 0; j < nen; j++){
                        for(int k = 0; k < nen; k++){
                            int ix = (iex + k + nx1_) % nx1_;
                            int iy = iey + j;
                            int iz = iez + i;
                            int id = iz * nx1_ * ny1_ + iy * nx1_ + ix;
                            ien.push_back(id);
                        }
                    }
                }
            }
        }
    }
    for(int iez = 0; iez < nez2; iez++){
        for(int iey = 0; iey <ney2; iey++){
            for(int iex = 0; iex < nex2; iex++){
                for(int i = 0; i < nen; i++){
                    for(int j = 0; j < nen; j++){
                        for(int k = 0; k < nen; k++){
                            int ix = iex + k;
                            int iy = iey + j;
                            int iz = iez + i;
                            int id = iz * nx2_ * ny2_ + iy * nx2_ + ix;   
                            ien.push_back(nn1+id);
                        }
                    }
                }
            }
        }
    }

    
    for(int ie = 0; ie <= ne; ie++){
        ien_offset.push_back(ie * 64);
    }
}

vector<double> IGA::bernstein_basis_function(double u)
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

vector<double>IGA::set_open_knot(int n)
{
    vector<double> knot;
	double int_knot = 0.0;
	for (int i = 0; i <= p; i++)
	{
		knot.push_back(int_knot);
	}
	int num_int_knot = n + p + 1 - 2 * (p + 1);
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

vector<double> IGA::bspline_basis_function(int n, double u,vector<double> knot)
{
	vector<double> N(n, 0.0);
	vector<double> Nn(n, 0.0);
	for (int i = 0; i < n; i++)
	{
		if (u >= knot.at(i) && u < knot.at(i + 1))
		{
			N.at(i) = 1.0;
		}
		else
		{
			N.at(i) = 0.0;
		}
	}
	if (u == knot.back())
	{
		N.at(n - 1) = 1.0;
	}
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double den1 = (knot.at(j + i + 1) - knot.at(j));
			double den2 = (knot.at(j + i + 2) - knot.at(j + 1));
			double turm1 = 0.0;
			double turm2 = 0.0;
			if (den1 != 0.0)
			{
				turm1 = (u - knot.at(j)) / den1 * N.at(j);
			}
			if (j < n - 1 && den2 != 0.0)
			{
				turm2 = (knot.at(j + i + 2) - u) / den2 * N.at(j + 1);
			}
			Nn.at(j) = turm1 + turm2;
		}
		for (int j = 0; j < n; j++)
		{
			N.at(j) = Nn.at(j);
		}
	}
	return N;
}

vector<double> IGA::bspline_curve(int nx)
{
	int np = 101;
	vector<double> C(2 * np, 0.0);
	vector<double> knot = IGA::set_open_knot(nx);
	double ks = knot.front();
	double ke = knot.back();
	for (int i = 0; i < np; i++)
	{
		double u = ks + (double)(ke - ks) * i / (np - 1);
		vector<double> N = IGA::bspline_basis_function(nx, u, knot);
		double cx = 0.0;
		double cy = 0.0;
		for (int j = 0; j < nx; j++)
		{
			cx += N.at(j) * cp.at(j);
			cy += N.at(j) * cp.at(nx + j);
		}
		C.at(i) = cx;
		C.at(np + i) = cy;
	}
	return C;
}

vector<double> IGA::bspline_surface(int nx, int ny)
{
	vector<double> C(2 * np * np, 0.0);
	vector<double> knot_x = IGA::set_open_knot(nx);
	vector<double> knot_y = IGA::set_open_knot(ny);
	double ksx = knot_x.front();
	double kex = knot_x.back();
	double ksy = knot_y.front();
	double key = knot_y.back();
	for (int i = 0; i < np; i++)
	{
		double uy = ksy + (key - ksy) * i / (np - 1);
		vector<double> Ny = bspline_basis_function(ny, uy, knot_y);
		for (int j = 0; j < np; j++)
		{
			double ux = ksx + (kex - ksx) * j / (np - 1);
			vector<double> Nx = bspline_basis_function(nx, ux, knot_x);
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

vector<double> IGA::bspline_volume(int nx, int ny, int nz)
{
	vector<double> C(3 * np * np * np, 0.0);
	vector<double> knot_x = IGA::set_open_knot(nx);
	vector<double> knot_y = IGA::set_open_knot(ny);
	vector<double> knot_z = IGA::set_open_knot(nz);
	double ksx = knot_x.front();
	double kex = knot_x.back();
	double ksy = knot_y.front();
	double key = knot_y.back();
	double ksz = knot_z.front();
	double kez = knot_z.back();
	for (int i = 0; i < np; i++)
	{
		double uz = ksz + (kez - ksz) * i / (np - 1);
		vector<double> Nz = bspline_basis_function(nz, uz, knot_z);
		for (int j = 0; j < np; j++)
		{
			double uy = ksy + (key - ksy) * j / (np - 1);
			vector<double> Ny = bspline_basis_function(ny, uy, knot_y);
			for (int k = 0; k < np; k++)
			{
				double ux = ksx + (kex - ksx) * k / (np - 1);
				vector<double> Nx = bspline_basis_function(nx, ux, knot_x);
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

vector<double> IGA::nurbs_curve(int nx)
{
	vector<double> C(2 * np, 0.0);
	vector<double> knot = IGA::set_open_knot(nx);
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
		vector<double> N = IGA::bspline_basis_function(nx, u, knot);
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

vector<double> IGA::nurbs_surface(int nx, int ny)
{
	vector<double> C(2 * np * np, 0.0);
	vector<double> knot_x = IGA::set_open_knot(nx);
	vector<double> knot_y = IGA::set_open_knot(ny);

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
		vector<double> Ny = IGA::bspline_basis_function(ny, uy, knot_y);
		for (int j = 0; j < np; j++)
		{
			double ux = ksx + (double)(kex - ksx) * j / (np - 1);
			vector<double> Nx = IGA::bspline_basis_function(nx, ux, knot_x);
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

vector<double> IGA::nurbs_volume(int nx, int ny, int nz)
{
	vector<double> C(3 * np * np * np, 0.0);
	vector<double> knot_x = IGA::set_open_knot(nx);
	vector<double> knot_y = IGA::set_open_knot(ny);
	vector<double> knot_z = IGA::set_open_knot(nz);

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
		vector<double> Nz = IGA::bspline_basis_function(nz, uz, knot_z);
		for (int j = 0; j < np; j++)
		{
		double uy = ksy + (double)(key - ksy) * j / (np - 1);
		vector<double> Ny = IGA::bspline_basis_function(ny, uy, knot_y);
			for(int k = 0; k < np; k++){
				double ux = ksx + (double)(kex - ksx) * k / (np - 1);
				vector<double> Nx = IGA::bspline_basis_function(nx, ux, knot_x);
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

vector<int> IGA::set_knotspan(vector<double> & knot)
{
	vector<int> knotspan;
	double tmp = knot.front();
	for (int i = 0; i < knot.size(); i++)
	{
		if (tmp != knot.at(i))
		{
			knotspan.push_back(i);
			tmp = knot.at(i);
		}
	}
	return knotspan;
}

vector<double> IGA::set_insert_knot(vector<double> & knot, vector<int> & knotspan)
{
	vector<double> insert_knot;
	for (int i = 0; i < knotspan.size(); i++)
	{
		int a = knotspan.at(i);
		double s = knot.at(a);
		for (int j = 0; j < p; j++)
		{
			if (knot.at(a + j) != s)
			{
				insert_knot.push_back(s);
			}
		}
	}
	return insert_knot;
}


vector<vector<double>> IGA::knot_insertion(int n){
	vector<double> knot = set_open_knot(n);
	vector<int> knotspan = set_knotspan(knot);
	vector<double> insert_knot = set_insert_knot(knot,knotspan);
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
    vector<vector<double>> C_ele(N-p,vector<double>(nen*nen,0.0));
    for(int ie = 0; ie < n-p; ie++){
        for(int i = 0; i < nen; i++){
            for(int j = 0; j < nen; j++){
                int row = ie + i;
                int col = ie * p + j;
                C_ele.at(ie).at(i*nen+j) = c.at(row * N + col);
            }
        }
    }

    return C_ele;
}

void IGA::periodic_tube(int nx){
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

    
    vector<double> w(nx*nen,0.0);
    for(int i = 0; i < w.size(); i++){
        w.at(i) = 1.0;
    }

    double pi = acos(-1.0);
    vector<double> tx;
    vector<double> ty;
    for(int i=0; i<nx; i++){
        double ang = 2.0*pi* static_cast<double>(i)/static_cast<double>(nx);
        double tmp_x = cos(ang);
        double tmp_y = sin(ang);
        tx.push_back(tmp_x);
        ty.push_back(tmp_y);
    }
    int nt = 10;
    double dt = (double)1.0/nt;
    double r_ave = 0.0;
    int count = 0;
    vector<double> p_x;
    vector<double> p_y;
    for(int i = 0; i < nx; i++){
        for(int it = 0; it < nt + 1; it++){
            double t = it * dt;
            vector<double> B = IGA::bernstein_basis_function(t);
            vector<double> R(nx*nen,0.0);
            for(int j = 0; j < nen; j++){
                for(int k = 0; k < nen; k++){
                    R.at(j) += C.at(j*nen+k)*B.at(k);
                }
            }
            double N = 0.0;
    
            for(int j = 0; j < nen; j++){
                N += w.at(j) * R.at(j);
            }

            double px =0.0;
            double py =0.0;
            for(int j = 0; j < nen; j++){
                int id = (i + j) % nx;
                px += tx.at(id) * R.at(j) / N;
                py += ty.at(id) * R.at(j) / N;
            }

            r_ave += sqrt(px*px + py*py); 
            count++;
        }
    }

    r_ave /= static_cast<double>(count);

    for(int i = 0; i < nn1; i++){
        cp.at(i) /= r_ave;
        cp.at(nn + i) /= r_ave;
    }
}

vector<double> IGA::periodic_circle_nurbs(int n, double u){
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

    vector<double> B = IGA::bernstein_basis_function(u);
    for(int i = 0; i < nen; i++){
        for(int j = 0; j < nen; j++){
            R.at(i) += C.at(i * nen + j) * B.at(j);
        }
    }
    return R;
}

vector<double> IGA::bezier_element(int n,double u,int ie){
    vector<double> R(p+1,0.0);
    vector<vector<double>> C = IGA::knot_insertion(n);
    vector<double> knot = IGA::set_open_knot(n);
    vector<int> knotspan = IGA::set_knotspan(knot);
    vector<double> insert_knot = IGA::set_insert_knot(knot, knotspan);
    int m = insert_knot.size();
    int N = n+m;

    vector<double> B = IGA::bernstein_basis_function(u);

    for(int i = 0; i < p + 1; i++){
            for(int j = 0; j < p + 1; j++){
                R.at(i) += C.at(ie).at(i*(p+1)+j) * B.at(j);
        }
    }
    return R;
}

vector<double> IGA::nurbs_basis(int ie,bool tube){
    int ien_offset = this->ien_offset.at(ie);
    int nex;
    int ney;
    int nez;

    int iex;
    int iey;
    int iez;

    int nx;
    int ny;
    int nz;

    if(tube){
        nx = nx1;
        ny = ny1;
        nz = nz1;
        nex = nx;
        ney = ny - p;
        nez = nz - p;
        iex = ie % nex;
        iey = (ie / nex) % ney;
        iez = ie / (nex * ney);
    }
    else{
        nx = nx2;
        ny = ny2;
        nz = nz2;
        nex = nx - p;
        ney = ny - p;
        nez = nz - p;
        iex = (ie-ne1) % nex;
        iey = ((ie-ne1) / nex) % ney;
        iez = (ie-ne1) / (nex * ney);
    }

    vector<double> nurbs;
    vector<double> w(nn,1.0);
    // for(int i = 0; i < nz; i++){
        // for(int j = 0; j < ny; j++){
            // for(int k = 0; k <nx; k++){
                // w.at(i*nx*ny+j*nx+k) = 1.0;
            // }
        // }
    // }

    for(int itz = 0; itz < np; itz++){
        double tz = (double)itz/(np-1); 
        vector<double> Rz = IGA::bezier_element(nz,tz,iez);
        for(int ity = 0; ity < np; ity++){
            double ty = (double)ity/(np-1);
            vector<double> Ry = IGA::bezier_element(ny,ty,iey);
            for(int itx = 0; itx < np; itx++){
                double tx = (double)itx/(np-1);
                vector<double> Rx;
                if(tube){
                     Rx = IGA::periodic_circle_nurbs(nx,tx);
                }
                else{
                     Rx = IGA::bezier_element(nx,tx,iex);
                }
                vector<double> R(nen*nen*nen, 0.0);
                for(int iz = 0; iz < nen; iz++){
                    for(int iy = 0; iy < nen; iy++){
                        for(int ix = 0; ix < nen; ix++){  
                            int id = iz*nen*nen + iy*nen + ix;
                            R.at(id) = Rx.at(ix) * Ry.at(iy) * Rz.at(iz);
                        }
                    }
                }
                double W = 0.0;
                for(int i = 0; i < nen; i++){
                    for(int j = 0; j < nen; j++){
                        for(int k = 0; k < nen; k++){
                            int id = i*nen*nen + j*nen + k;
                            int ien = this->ien.at(ien_offset+id);
                            W += w.at(ien) * R.at(id);
                        }
                    }
                }
                for(int i = 0; i < nen; i++){
                    for(int j = 0; j < nen; j++){
                        for(int k = 0; k < nen; k++){
                            int id = i*nen*nen+j*nen+k;
                            int ien = this->ien.at(ien_offset+id);
                            double n = w.at(ien) * R.at(id)/W;
                            nurbs.push_back(n);
                        }
                    }
                }
            }
        }
    }
    return nurbs;
}

vector<double> IGA::nurbs_iga(){
    #if 0 
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
            vector<double> N = IGA::periodic_circle_nurbs(nx, u);
    		double R = 0.0;
    		for (int j = 0; j < p+1; j++)
    		{
                int id = (ie + j + nx - 1) % nx;
    			R += w.at(id)* N.at(j);
    		}
    		double cx = 0.0;
    		double cy = 0.0;
    		for (int j = 0; j < p + 1; j++)
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
    int ney = ny -p;
    int ne = nex * ney;
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
	        	vector<double> Ny = IGA::bezier_element(ny,uy,iey);
	        	for (int j = 0; j < np; j++)
	        	{
                  double ux = (double)1.0/(np-1)*j;
	        		vector<double> Nx = IGA::periodic_circle_nurbs(nx,ux);
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
    vector<double> C(3 * ne * np * np * np, 0.0);
    IGA::periodic_tube(nx1);


    for(int ie = 0; ie < ne; ie++){
        vector<double> N;
        if(ie<ne1){
            bool tube = true;
            N = IGA::nurbs_basis(ie,tube);
        }
        else{
            bool tube = false;
            N = IGA::nurbs_basis(ie,tube);
        }
        int ien_offset_e = this->ien_offset.at(ie);
        for(int itz = 0; itz < np; itz++){
            for(int ity = 0; ity < np; ity++){
                for(int itx = 0; itx < np; itx++){
                    double cx = 0.0;
                    double cy = 0.0;
                    double cz = 0.0;
                    for(int iz = 0; iz < nen; iz++){
                        for(int iy = 0; iy < nen; iy++){
                            for(int ix = 0; ix < nen; ix++){
                                int id = itz*np*np*nen*nen*nen + ity*np*nen*nen*nen + itx*nen*nen*nen + iz*nen*nen + iy*nen + ix;
                                int ied = iz*nen*nen + iy*nen + ix;
                                int ien_e = this->ien.at(ien_offset_e+ied);
                                cx += N.at(id) * cp.at(ien_e);
                                cy += N.at(id) * cp.at(nn + ien_e);
                                cz += N.at(id) * cp.at(2 *nn + ien_e);
                            }
                        }
                    }
                    C.at(ie*np*np*np + itz*np*np + ity*np + itx) = cx;
                    C.at(ne*np*np*np + ie*np*np*np + itz*np*np + ity*np + itx) = cy;
                    C.at(2*ne*np*np*np + ie*np*np*np + itz*np*np + ity*np + itx) = cz;
                }
            }
        }
    }

    #endif
    return C;
}