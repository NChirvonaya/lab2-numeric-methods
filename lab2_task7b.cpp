#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

#include "lab1.h"

using namespace std;
vector<double> xx, a, b, c;

double f(double x)
{
  return pow(x, 3) + 3 * pow(x, 2) + 5;
}

double df(double x1, double x2)
{
  return f(x2) - f(x1);
}

double ddf(double x1, double x2, double x3)
{
  return df(x2, x3) - df(x1, x2);
}

double p(int k, double x)
{
  return a[k] + b[k] * (x - xx[k]) * c[k - 1] * pow(x - xx[k], 2);
}

int main()
{
  
  freopen("out.plot", "w", stdout);

  int n;
  cin >> n;
  double l, r;
  cin >> l >> r;
  vector<double> x(n + 1);
  xx.resize(n + 2);
  a.resize(n + 2);
  b.resize(n + 2); 
  c.resize(n + 1);
  double dd(1.0*(r - l) / n);
  x[0] = l;
  for (int i(1); i <= n; ++i) {
    x[i] = x[i-1] + i*dd;
    xx[i] = (i == 0 ? x[i] : (x[i - 1] + x[i]) / 2);
  }
  xx[n + 1] = x[n];
  
  vector<double> h(n + 1);
  for (int i(1); i <= n; ++i) {
    h[i] = x[i] - x[i - 1];
  }
  
  double A(0), B(0);

  vector<vector<double> > mat(n + 1, vector<double>(n + 1, 0));
  vector<double> fr(n + 2);
  mat[0][0] = 2;
  fr[0] = A;
  mat[1][1] = 3;
  mat[1][2] = h[2]/(h[1]+h[2]);
  fr[1] = 4 * ddf(x[0], x[1], x[2]) - A * h[1] / (2 * h[1] + 2 * h[2]);
  for (int i(2); i < n - 1; ++i) {
    mat[i][i-1] = h[i] / (h[i] + h[i + 1]);
    mat[i][i] = 3;
    mat[i][i + 1] = h[i + 1] / (h[i] + h[i + 1]);
    fr[i] = 4 * ddf(x[i - 1], x[i], x[i + 1]);
  }
  mat[n - 1][n - 2] = h[n - 1] / (h[n - 1] + h[n]);
  mat[n - 1][n-1] = 3;
  fr[n - 1] = 4 * ddf(x[n - 2], x[n - 1], x[n]) - B*h[n] / (2 * h[n - 1] + 2 * h[n]);
  mat[n][n] = 2;
  fr[n] = B;

  SLAE sys(mat, fr);
  sys.rotationMethod();
  vector<double> c = sys.solution;

  for (int i(1); i <= n; ++i) {
    b[i] = h[i] * (c[i - 1] - c[i])/4 + df(x[i - 1], x[i]);
  }
  b[n + 1] = b[n] + B*h[n] / 2;

  for (int i(1); i <= n; ++i) {
    a[i] = h[i] * (b[i] - h[i] * c[i - 1] / 2) / 2 + f(x[i - 1]);
  }
  a[n + 1] = f(x[n]);

  for (int i(0); i <= n; ++i) {
    printf("%5.3lf %5.3lf ", x[i], f(x[i]));
    printf("%5.3lf\n", p(i + 1, x[i]));
  }
  return 0;
}