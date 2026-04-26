#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> bernstein(int p, double u)
{
  vector<double> B(p + 1, 0.0);
  vector<double> Bn(p + 1, 0.0);
  for (int i = 0; i < p + 1; i++)
  {
    B.at(i) = 1.0;
  }
  for (int i = 0; i < p ; i++)
  {
    Bn.at(0) = (1 - u) * B.at(0);
    for (int j = 1; j < p + 1; j++)
    {
      Bn.at(j) = (1 - u) * B.at(j) + u * B.at(j - 1);
    }
    for (int j = 0; j < p + 1; j++)
    {
      B.at(j) = Bn.at(j);
    }
  }
  return B;
}

int main()
{
  int p = 3;
  vector<double> cp(2 * (p + 1), 0.0);
  cp.at(0) = 1.0;
  cp.at(1) = 2.5;
  cp.at(2) = 3.4;
  cp.at(3) = 4.8;
  cp.at(4) = 1.3;
  cp.at(5) = 2.6;
  cp.at(6) = 3.5;
  cp.at(7) = 2.6;

  vector<double> B;

  for (int i = 0; i < 100; i++)
  {
    double u = i / 100.0;
    B = bernstein(p, u);

    vector<double> C(2, 0.0);

    for (int a = 0; a < p + 1; a++)
    {
      C.at(0) += B.at(a) * cp.at(a);
      C.at(1) += B.at(a) * cp.at(p + 1 + a);
    }
    cout << C.at(0) << ", " << C.at(1) << endl;
  }
  return 0;
}