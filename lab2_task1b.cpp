#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

#include "lab1.h"

using namespace std;
vector<double> c;
double m;

double f(double x)
{
  return log(pow(x, 3)) + 3 * pow(x, 2) + 5;
}

double p(double x)
{
  double res(0.0);
  for (int i(0); i <= m; ++i) {
    res += c[i] * pow(x, i);
  }
  return res;
}

int main()
{

  int n;
  cin >> n;
  double l, r;
  cin >> l >> r;
  vector<double> x(n + 1);
  double dd(1.0*(r - l) / n);
  x[0] = l;
  for (int i(1); i <= n; ++i) {
    x[i] = x[i - 1] + i*dd;
  }

  cin >> m;

  c.resize(m + 1);

  vector<vector<double> > mat(m + 1, vector<double>(m + 1, 0.0));
  vector<double> fr(m + 1, 0.0);
  for (int i(0); i <= m; ++i) {
    for (int j(0); j <= m; ++j) {
      for (int k(0); k <= n; ++k) {
        mat[i][j] += pow(x[k], i + j);
      }
      mat[i][j] /= (n + 1);
    }
    for (int k(0); k <= n; ++k) {
      fr[i] += pow(x[k], i)*f(x[k]);
    }
    fr[i] /= (n + 1);
  }

  SLAE sys(mat, fr);
  sys.rotationMethod();
  c = sys.solution;

  for (int i(0); i <= n; ++i) {
    printf("%8.3lf %8.3lf ", x[i], f(x[i]));
    printf("%8.3lf\n", p(x[i]));
  }
  return 0;
}