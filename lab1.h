#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

const double eps = 1e-3;

std::vector<double> operator*(const std::vector<double>& v, double k);

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2);

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);

std::vector<double> operator/(const std::vector<double>& v1, const std::vector<double>& v2);

std::vector<double> operator*(const std::vector<std::vector<double> >& A, const std::vector<double>& v);

std::vector<std::vector<double> > operator*(const std::vector<std::vector<double> >& v1, const std::vector<std::vector<double> >& v2);

std::vector<std::vector<double> > operator+(const std::vector<std::vector<double> >& v1, const std::vector<std::vector<double> >& v2);

std::vector<std::vector<double> > operator-(const std::vector<std::vector<double> >& v1, const std::vector<std::vector<double> >& v2);

std::vector<std::vector<double> > operator*(double k, const std::vector<std::vector<double> >& v);

std::vector<double> normalize(std::vector<double> v);

double norm2(const std::vector<double>& v);

double average(const std::vector<double>& v);

double det(const std::vector<std::vector<double> >& A);

double scal(const std::vector<double>& v1, const std::vector<double>& v2);

bool checkConvergence(const std::vector<double>& v1, const std::vector<double>& v2);

class SLAE
{
public:
	std::vector<std::vector<double> > A;
	std::vector<std::vector<double> > L;
	std::vector<std::vector<double> > U;
	std::vector<double> B;
	std::vector<double> solution;
	std::vector<double> eigenValues;
	std::vector<std::vector<double> > eigenVectors;

	int varcount;
	int eqcount;
	SLAE(const std::vector<std::vector<double> >& _A, const std::vector<double>& _B);

	SLAE(const std::vector<std::vector<double> >& _A);
	
	void reverseRun(int k);

	void rotationMethod();

	void LUdecomposition();

	void LUsolution();

	std::vector<std::vector<double> > inverseMatrix();

	void invMatMethod();

	void sqRootsMethod();

	void eigenValuesLU();

	void eigenVector(int i);

	std::vector<double> eigenVector(double lambda);

	void conjugateGradientMethod();

	double maxEigenValue();

	double minEigenValue();

};


