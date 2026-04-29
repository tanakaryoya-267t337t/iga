#include <iostream>
#include <vector>
#include <cmath>
#include "output.h"
#include "iga.h"

using namespace std;

vector<double> bernstein(int p, double u)
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

vector<double> bspline(int p, int a, double u, vector<double> &k)
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

void knot_insert(int p, int a, vector<double> &knot, vector<double> &cp, double u)
{
	int n = knot.size();
	int ip = 0;
	vector<double> new_knot(n+1,0.0);
	for (int k = 0; k < n-1; k++)
	{
		if (u >= knot.at(k) && u < knot.at(k + 1))
		{
			ip = k;
			new_knot.at(ip + 1) = u;
			int N = n - (ip + 1);
			for (int j = 0; j < n+1; j++)
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

	for (int i = 0; i < a + 1; i++)
	{
		if (i <= ip - p)
		{
			CP.at(i) = cp.at(i);
			CP.at(a + 1 + i) = cp.at(a + i);
		}
		else if (i > ip - p && i <= ip)
		{
			double alpha = (double)(u - knot.at(i)) / (knot.at(i + p) - knot.at(i));
			if((u-knot.at(i)) == 0){
				alpha = 0.0;
			}
			CP.at(i) = alpha * cp.at(i) + (1.0 - alpha) * cp.at(i - 1);
			CP.at(a + 1 + i) = alpha * cp.at(a + i) + (1.0 - alpha) * cp.at(a + i - 1);
		}
		else
		{
			CP.at(i) = cp.at(i - 1);
			CP.at(a + 1 + i) = cp.at(a + i - 1);
		}
	}
	knot = new_knot;
	cp = CP;
}
// void knot_insertion(
// 		vector<double> &knotvector,
// 		vector<double> &insert_knot,
// 		vector<vector<double>> &Cm,
// 		int ncp,
// 		int p)
// {
// 	double ins = insert_knot.front();
// 	int k;
// 	for (int i = 0; i < knotvector.size() - 1; i++)
// 	{
// 		if (knotvector.at(i) <= ins && ins < knotvector.at(i + 1))
// 		{
// 			k = i; // place to insert knot
// 			break;
// 		}
// 	}
// 	int m = ncp + 1;

// 	Cm.resize(ncp);
// 	for (int i = 0; i < Cm.size(); i++)
// 	{
// 		Cm.at(i).resize(m);
// 		for (int j = 0; j < Cm.at(i).size(); j++)
// 		{
// 			Cm.at(i).at(j) = 0.0;
// 		}
// 	}

// 	for (int i = 0; i < m; i++)
// 	{
// 		double alpha;
// 		if (i <= k - p)
// 		{
// 			alpha = 1.0;
// 		}
// 		else if (i <= k)
// 		{
// 			alpha = (ins - knotvector.at(i)) / (knotvector.at(i + p) - knotvector.at(i));
// 		}
// 		else
// 		{
// 			alpha = 0.0;
// 		}

// 		if (i < ncp)
// 		{
// 			Cm.at(i).at(i) = alpha;
// 		}
// 		if (i > 0)
// 		{
// 			Cm.at(i - 1).at(i) = 1.0 - alpha;
// 		}
// 	}

// 	knotvector.insert(knotvector.begin() + k + 1, ins); // knot(k), newknot, knot(k+1)
// 	insert_knot.erase(insert_knot.begin());
// }