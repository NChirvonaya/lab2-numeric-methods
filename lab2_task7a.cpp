#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;

const double eps = 1e-6;

double f(double x)
{
  return pow(x, 3) - 0.2 * pow(x, 2) - 0.2 * x - 1.2;
}

double phi(double x)
{
  return pow(0.2 * pow(x, 2) + 0.2 * x + 1.2, 1.0/3.0);
}

int main()
{
  freopen("out7.plot", "w", stdout);
/*
  cout << "y(x) = x^3 - 0.2x^2 - 0.2x - 1.2\n" << "y(x) = 0\n";

  cout << "Chord method: \n";*/
  double a(1), b(1.5);
  double q(0.170849);
  double x0(a);
  double x(x0 - f(x0) / (f(b) - f(x0))* (b - x0));
  while (fabs(x - x0) > eps) {
    x0 = x;
    x = x0 - f(x0) / (f(b) - f(x0))* (b - x0);
  }/*
  cout << "x0 = " << x << "\n";
  cout << "y(x0) = " << fixed << setprecision(5) << f(x) << "\n";


  cout << "Simple iterations method: \n";*/
  x0 = a;
  x = phi(x0);
  while (fabs(x - x0) > eps) {
    printf("%8.3lf %8.3lf\n", x, 0.0);
    x0 = x;
    x = phi(x0);
  }/*
  cout << "x0 = " << x << "\n";
  cout << "y(x0) = " << fixed << setprecision(5) << f(x) << "\n";

  cout << "Vegstein method: \n";*/
  x0 = a;
  double x1(phi(x0));
  double nx0(x0);
  double nx1(x1);
  double x2(phi(nx1));
  double nx2;
  while (fabs(x2 - nx1) > eps*(1.0 - q)/q) {
    nx2 = (x2 * nx0 - x1 * nx1) / (x2 + nx0 - x1 - nx1);
    nx0 = x1;
    x1 = x2;
    nx1 = nx2;
  }/*
  cout << "x0 = " << x1 << "\n";
  cout << "y(x0) = " << fixed << setprecision(5) << f(x) << "\n";*/

  return 0;
}