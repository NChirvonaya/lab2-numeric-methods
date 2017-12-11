#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;

const double eps = 1e-6;

double f(double x, double y)
{
  return pow(x, 2) - 2 * pow(y, 2) - y * x + 2 * x - y + 1;
}

double g(double x, double y)
{
  return 2 * pow(x, 2) - pow(y, 2) + y * x + 3 * y - 5;
}

double df_dx(double x, double y)
{
  return 2 * x - y + 2;
}

double dg_dx(double x, double y)
{
  return 4 * x + y;
}

double df_dy(double x, double y)
{
  return -4 * y - x - 1;
}

double dg_dy(double x, double y)
{
  return -2 * y + x + 3;
}
int main()
{
  cout << "f(x, y) = x^2 - 2y^2 - xy + 2x - y + 1\n";
  cout << "g(x, y) = 2x^2 - y^2 - yx + 3y - 5\n";

  cout << "Brown method:\n";

  double xk(0), x_k(0), qk(1), pk(1), yk(0);

  while (max(abs(pk), abs(qk)) > eps) {
    x_k = xk - f(xk, yk) / df_dx(xk, yk);
    qk = g(x_k, yk) * df_dx(xk, yk) / (df_dx(xk, yk)*dg_dy(xk, yk) - df_dy(xk, yk)*dg_dx(xk, yk));
    pk = (f(xk, yk) - qk*df_dy(xk, yk)) / df_dx(xk, yk);
    xk -= pk;
    yk -= qk;
  }

  cout << fixed << setprecision(5) << "x0 = " << xk << "\n" << "y0 = " << yk << "\n";
  cout << "f(x0,y0) = " << fixed << setprecision(5) << f(xk, yk) << "\n";
  cout << "g(x0,y0) = " << fixed << setprecision(5) << g(xk, yk) << "\n";

  return 0;
}