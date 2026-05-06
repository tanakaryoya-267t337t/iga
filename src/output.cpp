#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include "output.h"

using namespace std;

void Output(int np,int nn, int ne, vector<double> &x, string &filename)
{
	ofstream file("result/vtk/" + filename);
	file << "# vtk DataFile Version 2.0" << endl;
	file << "Title Data" << endl;
	file << "ASCII" << endl;
	file << "DATASET UNSTRUCTURED_GRID" << endl;
	file << "POINTS " << nn << " double" << endl;
	#if 1 
	for (int i = 0; i < nn; i++)
	{
		file << x.at(i) << " " << x.at(nn + i) << " 0" << endl;
	}
	#if 1 
	file << "CELLS " << ne << " " << 3 * ne << endl;
	for (int i = 0; i < ne; i++)
	{
		file << "2 " << i << " " << i + 1 << endl;
	}
	file << "CELL_TYPES " << ne << endl;
	for (int i = 0; i < ne; i++)
	{
		file << "3" << endl;
	}
	#else 
  file << "CELLS " << ne << " " << 5 * ne << endl;
  for (int i = 0; i < np - 1; i++)
  {
      for(int j = 0; j < np - 1; j++)
      {
          int id = i * np + j;
          file << "4 " << id << " " << id + 1 << " " << id + np + 1 << " " << id + np << endl;
      }
  }
  file << "CELL_TYPES " << ne << endl;
  for (int i = 0; i < ne; i++)
  {
      file << "9" << endl;
  }
	#endif
	#else
	for (int i = 0; i < nn; i++)
	{
		file << x.at(i) << " " << x.at(nn + i) << " " << x.at(2*nn + i) << endl;
	}
  file << "CELLS " << ne << " " << 9 * ne << endl;
  for (int i = 0; i < np - 1; i++)
  {
      for(int j = 0; j < np - 1; j++)
      {
				for(int k = 0; k < np - 1; k++){
          int id = i * np * np + j * np + k;
          file << "8 " << id << " " << id + 1 << " " << id + np + 1 << " " << id + np << " " << id + np * np << " " << id + np * np + 1 << " " << id + np * np + np + 1 << " " << id + np * np + np << endl;
				}
      }
  }
  file << "CELL_TYPES " << ne << endl;
  for (int i = 0; i < ne; i++)
  {
      file << "12" << endl;
  }
	#endif
	file.close();
}

void output_gnuplot(int n, int nn, vector<vector<double>> &B, string &gnuplotname)
{
	FILE *gp = popen("gnuplot", "w");
	fprintf(gp, "set terminal svg size 800,600\n");
	fprintf(gp, "set output 'result/gnuplot/%s'\n", gnuplotname.c_str());
	fprintf(gp, "set object 1 rect from screen 0,0 to screen 1,1 behind fc rgb 'white' fillstyle solid 1.0\n");fprintf(gp, "set xtics (0, 1, 2, 3, 4)\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "set xlabel 'ξ'\n");
	fprintf(gp, "set ylabel 'N(ξ)'\n");
	fprintf(gp, "set xtics nomirror\n");
	fprintf(gp, "set ytics nomirror\n");

	fprintf(gp, "plot ");
	for (int a = 0; a < n; a++)
	{
		fprintf(gp, "'-' w l title 'N%d'", a);
		if (a != n - 1)
			fprintf(gp, ", ");
	}
	fprintf(gp, "\n");

	for (int a = 0; a < n; a++)
	{
		for (int i = 0; i < nn; i++)
		{
			double u = 4.0 * i / (double)(nn - 1);
			fprintf(gp, "%f %f\n", u, B.at(a).at(i));
		}
		fprintf(gp, "e\n");
	}
	pclose(gp);
}