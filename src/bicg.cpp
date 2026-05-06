#include <cmath>
#include <iostream>
#include <vector>
#include "bicg.h"

using namespace std;

// 内積
double dot(vector<double>& a, vector<double>& b) {
  double result = 0.0;
  for (int i = 0; i < a.size(); i++) {
    result += a[i] * b[i];
  }
  return result;
}

// 行列の積
double product(vector<double>& A, vector<double>& b, int i) {
  int n         = b.size();
  double result = 0.0;
  for (int j = 0; j < n; j++) {
    result += A[i * n + j] * b[j];
  }
  return result;
}

// 2ノルム
double norm(vector<double>& a) {
  return sqrt(dot(a, a));
}

vector<double> matT(vector<double>& A, int col, int row) {
  vector<double> AT(row * col, 0.0);
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      AT[j * row + i] = A[i * col + j];
    }
  }
  return AT;
}

vector<double> bicg(vector<double>& A_, vector<double>& b_, int N, bool bc) {
  vector<double> u(N, 0.0);
  vector<double> A;
  vector<double> b;
  int n;
  vector<double> AT;
  vector<double> x;
  if (bc == true) {
    const int Nx = 5;
    const int Ny = 4;
    const int ny = Ny - 2;
    const int n  = Nx * ny;

    for (int i = 0; i < Nx; i++) {
      u[i]                 = 0.0;
      u[Nx * (Ny - 1) + i] = 1.0;
    }

    x.assign(n, 0.0);

    vector<double> A(n * n, 0.0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i * n + j] = A_[(i + 1) * N + (j + 1)];
      }
    }
    vector<double> AT(n * n, 0.0);

    AT = matT(A, n, n);

    vector<double> b(n, 0.0);
    for (int i = 0; i < n; i++) {
      b[i] = b_[i + 1];
    }
  } else {
    n = N;
    x.assign(n, 0.0);
    A  = A_;
    AT = matT(A, n, n);
    b  = b_;
  }
  vector<double> r(n, 0.0);
  vector<double> r_tld(n, 0.0);
  vector<double> p(n, 0.0);
  vector<double> p_tld(n, 0.0);

  for (int i = 0; i < n; i++) {
    r[i]     = b[i] - product(A, x, i);
    p[i]     = r[i];
    r_tld[i] = r[i];
    p_tld[i] = r_tld[i];
  }

  const int iter_max = 1000;
  const double Tol   = 1.0e-20;

  for (int i = 0; i < iter_max; i++) {
    double alpha = 0.0;
    double beta  = 0.0;
    vector<double> A_p(n, 0.0);
    vector<double> p_tld_A(n, 0.0);
    vector<double> AT_p_tld(n, 0.0);
    for (int j = 0; j < n; j++) {
      A_p[j]      = product(A, p, j);
      p_tld_A[j]  = product(A, p_tld, j);
      AT_p_tld[j] = product(AT, p_tld, j);
    }
    double denominator = dot(p_tld, A_p);
    alpha              = dot(r_tld, r) / denominator;

    vector<double> r_n(n, 0.0);
    vector<double> r_tld_n(n, 0.0);
    vector<double> x_n(n, 0.0);

    for (int j = 0; j < n; j++) {
      x_n[j]     = x[j] + alpha * p[j];
      r_n[j]     = r[j] - alpha * A_p[j];
      r_tld_n[j] = r_tld[j] - alpha * AT_p_tld[j];
    }

    double denominator_beta = dot(r_tld, r);

    beta = dot(r_tld_n, r_n) / denominator_beta;

    vector<double> p_n(n, 0.0);
    vector<double> p_tld_n(n, 0.0);

    for (int j = 0; j < n; j++) {
      p_n[j]     = r_n[j] + beta * p[j];
      p_tld_n[j] = r_tld_n[j] + beta * p_tld[j];
    }

    for (int j = 0; j < n; j++) {
      x[j]     = x_n[j];
      r[j]     = r_n[j];
      r_tld[j] = r_tld_n[j];
      p[j]     = p_n[j];
      p_tld[j] = p_tld_n[j];
    }

    double error = norm(r);

    if (error < Tol) {
      break;
    }
  }

  for (int i = 0; i < n; i++) {
    u[i] = x[i];
  }
  return u;
}