#include "lab1.h"

std::vector<double> operator*(const std::vector<double>& v, double k)
{
	std::vector<double> res;
	int n(v.size());
	for (int i(0); i < n; ++i) {
		res.push_back(v[i] * k);
	}
	return res;
}

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2)
{
	std::vector<double> res;
	int n1(v1.size());
	int n2(v2.size());
	try {
		throw n1 != n2;
	}
	catch (bool f) {
		if (f) {
			std::cout << "You can only sum equal-seized vectors!" << std::endl;
		}
	}

	for (int i(0); i < n1; ++i) {
		res.push_back(v1[i] + v2[i]);
	}
	return res;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
	std::vector<double> res;
	int n1(v1.size());
	int n2(v2.size());
	try {
		throw n1 != n2;
	}
	catch (bool f) {
		if (f) {
			std::cout << "You can only diff equal-seized vectors!" << std::endl;
		}
	}

	for (int i(0); i < n1; ++i) {
		res.push_back(v1[i] - v2[i]);
	}
	return res;
}

std::vector<double> operator/(const std::vector<double>& v1, const std::vector<double>& v2)
{
	std::vector<double> res;
	int n1(v1.size());
	int n2(v2.size());
	try {
		throw n1 != n2;
	}
	catch (bool f) {
		if (f) {
			std::cout << "You can only diff equal-seized vectors!" << std::endl;
		}
	}

	for (int i(0); i < n1; ++i) {
		res.push_back(v1[i] / v2[i]);
	}
	return res;
}

std::vector<double> operator*(const std::vector<std::vector<double> >& A, const std::vector<double>& v)
{
	try {
		throw A[0].size() != v.size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "Wrong matrix or vector size!" << std::endl;
		}
	}

	int n(A.size());
	int l(v.size());
	std::vector<double> res(n);

	for (int i(0); i < res.size(); ++i) {
		res[i] = 0;
		for (int j(0); j < l; ++j) {
			res[i] += A[i][j] * v[j];
		}
	}

	return res;
}

std::vector<std::vector<double> > operator*(const std::vector<std::vector<double> >& v1, const std::vector<std::vector<double> >& v2)
{
	try {
		throw v1[0].size() != v2.size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "Wrong matrixes size!" << std::endl;
		}
	}
	std::vector<std::vector<double> > res;
	int n(v1.size());
	int m(v2[0].size());
	int l(v2.size());
	res.resize(n);
	for (int i(0); i < res.size(); ++i) {
		res[i].resize(m);
	}

	for (int i(0); i < n; ++i) {
		for (int j(0); j < m; ++j) {
			res[i][j] = 0;
			for (int k(0); k < l; ++k) {
				res[i][j] += v1[i][k] * v2[k][j];
			}
		}
	}

	return res;
}

std::vector<std::vector<double> > operator+(const std::vector<std::vector<double> >& v1, const std::vector<std::vector<double> >& v2)
{
	try {
		throw v1.size() != v2.size() || v1[0].size() != v2[0].size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "Wrong matrixes size!" << std::endl;
		}
	}
	std::vector<std::vector<double> > res(v1.size());

	for (int i(0); i < v1.size(); ++i) {
		for (int j(0); j < v1[i].size(); ++j) {
			res[i].push_back(v1[i][j] + v2[i][j]);
		}
	}

	return res;
}

std::vector<std::vector<double> > operator-(const std::vector<std::vector<double> >& v1, const std::vector<std::vector<double> >& v2)
{
	try {
		throw v1.size() != v2.size() || v1[0].size() != v2[0].size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "Wrong matrixes size!" << std::endl;
		}
	}
	std::vector<std::vector<double> > res(v1.size());

	for (int i(0); i < v1.size(); ++i) {
		for (int j(0); j < v1[i].size(); ++j) {
			res[i].push_back(v1[i][j] - v2[i][j]);
		}
	}

	return res;
}

std::vector<std::vector<double> > operator*(double k, const std::vector<std::vector<double> >& v)
{
	std::vector<std::vector<double> > res(v.size());

	for (int i(0); i < v.size(); ++i) {
		for (int j(0); j < v[i].size(); ++j) {
			res[i].push_back(v[i][j] * k);
		}
	}

	return res;
}

std::vector<double> normalize(std::vector<double> v)
{
	double mod(0.0);
	for (int i(0); i < v.size(); ++i) {
		mod += v[i] * v[i];
	}

	mod = sqrt(mod);

	for (int i(0); i < v.size(); ++i) {
		v[i] /= mod;
	}

	return v;
}

double norm2(const std::vector<double>& v)
{
	double res(0.0);

	for (int i(0); i < v.size(); ++i) {
		res = std::max(res, std::abs(v[i]));
	}

	return res;
}

double average(const std::vector<double>& v)
{
	double res(0.0);

	for (int i(0); i < v.size(); ++i) {
		res += v[i];
	}

	res /= v.size();

	return res;
}

double det(const std::vector<std::vector<double> >& A)
{
	try {
		throw A[0].size() != A.size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "You can find determinant only for square matrix!" << std::endl;
		}
	}

	if (A.size() == 1)
		return A[0][0];

	double ans(0.0);

	int z = 1;
	for (int k(0); k < A[0].size(); ++k) {
		std::vector<std::vector<double> > a;
		a.resize(A.size() - 1);
		for (int i(0); i < a.size(); ++i) {
			a[i].resize(A[0].size() - 1);
		}
		for (int i(1); i < A.size(); ++i) {
			int jj(0);
			for (int j(0); j < A[i].size(); ++j) {
				if (j != k) {
					a[i - 1][jj] = A[i][j];
					jj++;
				}
			}
		}
		ans += z * det(a) * A[0][k];
		z *= -1;
	}
	return ans;
}

double scal(const std::vector<double>& v1, const std::vector<double>& v2)
{
	try {
		throw v1.size() != v2.size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "Wrong vectors size!" << std::endl;
		}
	}

	double res(0.0);

	for (int i(0); i < v1.size(); ++i) {
		res += v1[i] * v2[i];
	}

	return res;
}

bool checkConvergence(const std::vector<double>& v1, const std::vector<double>& v2)
{
	try {
		throw v1.size() != v2.size();
	}
	catch (bool f) {
		if (f) {
			std::cout << "Wrong vectors size!" << std::endl;
		}
	}

	int i(0);
	while (i < v1.size() && abs(v2[i] - v1[i]) < eps) {
		i++;
	}

	return i < v1.size();
}

SLAE::SLAE(const std::vector<std::vector<double> >& _A, const std::vector<double>& _B) :A(_A), B(_B), varcount(_A[0].size()), eqcount(_A.size()) {
	solution.resize(varcount);
}

SLAE::SLAE(const std::vector<std::vector<double> >& _A) :A(_A), varcount(_A[0].size()), eqcount(_A.size()) {
	solution.resize(varcount);
	B.assign(eqcount, 0);
}

void SLAE::reverseRun(int k)
{
	for (int i(k - 1); i >= 0; --i) {
		solution[i] = B[i];
		for (int j(i + 1); j < varcount; ++j) {
			solution[i] -= solution[j] * A[i][j];
		}
		solution[i] /= A[i][i];
	}
}

void SLAE::rotationMethod()
{
	for (int k(0); k < eqcount - 1; ++k) {
		for (int n(k + 1); n < eqcount; ++n) {
			double c(A[k][k] / sqrt(A[k][k] * A[k][k] + A[n][k] * A[n][k]));
			double s(A[n][k] / sqrt(A[k][k] * A[k][k] + A[n][k] * A[n][k]));
			std::vector<double> newAk;
			double newBk;
			std::vector<double> newAn;
			double newBn;
			newAk = (A[k] * c) + (A[n] * s);
			newBk = B[k] * c + B[n] * s;
			newAn = (A[k] * (-s)) + (A[n] * c);
			newBn = B[k] * (-s) + B[n] * c;
			A[k] = newAk;
			B[k] = newBk;
			A[n] = newAn;
			B[n] = newBn;
		}
	}
	reverseRun(eqcount);
}

void SLAE::LUdecomposition()
{
	L.resize(eqcount);
	U.resize(eqcount);
	for (int i(0); i < eqcount; ++i) {
		L[i].resize(varcount);
		U[i].resize(varcount);
	}

	for (int i(0); i < eqcount; ++i) {
		for (int j(0); j < varcount; ++j) {
			if (i <= j) {
				U[i][j] = A[i][j];
				for (int t(0); t <= i; ++t) {
					U[i][j] -= L[i][t] * U[t][j];
				}
			}
			else {
				L[i][j] = A[i][j];
				for (int t(0); t < j; ++t) {
					L[i][j] -= L[i][t] * U[t][j];
				}
				L[i][j] /= U[j][j];
			}
		}
		L[i][i] = 1;
	}
}

void SLAE::LUsolution()
{
	LUdecomposition();
	std::vector<double> y(varcount);
	for (int j(0); j < varcount; ++j) {
		y[j] = B[j];
		for (int i(0); i < j; ++i) {
			y[j] -= L[j][i] * y[i];
		}
	}

	for (int j(eqcount - 1); j >= 0; --j) {
		solution[j] = y[j];
		for (int i(j + 1); i < varcount; ++i) {
			solution[j] -= U[j][i] * solution[i];
		}
		solution[j] /= U[j][j];
	}
}

std::vector<std::vector<double> > SLAE::inverseMatrix()
{
	LUdecomposition();
	std::vector<std::vector<double> > inv;
	inv.resize(eqcount);
	for (int i(0); i < eqcount; ++i) {
		inv[i].resize(varcount);
	}

	std::vector<double> b;
	for (int i(0); i < eqcount; ++i) {
		b.clear();
		for (int j(0); j < eqcount; ++j) {
			if (j == i) {
				b.push_back(1);
			}
			else {
				b.push_back(0);
			}
		}
		SLAE s(A, b);
		s.LUsolution();
		for (int j(0); j < eqcount; ++j) {
			inv[j][i] = s.solution[j];
		}
	}
	return inv;
}

void SLAE::invMatMethod()
{
	std::vector < std::vector<double> > invA(inverseMatrix());
	solution = invA * B;
}

void SLAE::sqRootsMethod()
{
	std::vector<std::vector<double> > u;
	u.resize(eqcount);
	for (int i(0); i < u.size(); ++i) {
		for (int j(0); j < varcount; ++j) {
			u[i].push_back(0);
		}
	}

	for (int i(0); i < u.size(); ++i) {
		u[i][i] = A[i][i];
		for (int k(0); k < i; ++k) {
			u[i][i] -= (u[k][i] * u[k][i]);
		}
		u[i][i] = sqrt(u[i][i]);

		for (int j(i + 1); j < varcount; ++j) {
			u[i][j] = A[i][j];
			for (int k(0); k < i; ++k) {
				u[i][j] -= (u[k][i] * u[k][j]);
			}
			u[i][j] /= u[i][i];
		}
	}

	std::vector<double> y(varcount);
	for (int i(0); i < y.size(); ++i) {
		y[i] = B[i];
		for (int k(0); k < i; ++k) {
			y[i] -= u[k][i] * y[k];
		}
		y[i] /= u[i][i];
	}

	for (int i(varcount - 1); i >= 0; --i) {
		solution[i] = y[i];
		for (int k(i + 1); k < varcount; ++k) {
			solution[i] -= u[i][k] * solution[k];
		}
		solution[i] /= u[i][i];
	}
}

void SLAE::eigenValuesLU()
{
	SLAE a(A);
	eigenValues.assign(eqcount, 0);
	eigenVectors.resize(eqcount);

	std::vector<std::vector<double> > b(eqcount);
	int k(0);
	while (k < 20) {
		a.LUdecomposition();
		SLAE newA(a.U * a.L);
		for (int i(0); i < eqcount; ++i) {
			b[i].resize(varcount);
			for (int j(0); j < varcount; ++j) {
				b[i][j] = A[i][j];
				if (i == j) {
					b[i][j] -= newA.A[i][i];
				}
			}
			k++;
		}
		a = newA;
	}
	for (int i(0); i < eqcount; ++i) {
		eigenValues[i] = a.A[i][i];
		eigenVector(i);
	}
}

void SLAE::eigenVector(int k)
{
	try {
		throw k < 0 || k >= eqcount;
	}
	catch (bool f) {
		if (f) {
			std::cout << "Invalid argument!" << std::endl;
		}
	}
	std::vector<std::vector<double> > a(A);
	for (int i(0); i < eqcount; ++i) {
		a[i][i] -= eigenValues[k];
	}
	std::vector<double > b;
	for (int i(0); i < eqcount; ++i) {
		b.push_back(0);
	}

	SLAE S(a, b);
	S.solution[eqcount - 1] = 1;
	S.reverseRun(eqcount - 1);
	eigenVectors[k] = normalize(S.solution);
}

std::vector<double> SLAE::eigenVector(double lambda)
{
	std::vector<std::vector<double> > a(A);
	for (int i(0); i < eqcount; ++i) {
		a[i][i] -= lambda;
	}
	std::vector<double > b;
	for (int i(0); i < eqcount; ++i) {
		b.push_back(0);
	}

	SLAE S(a, b);
	S.solution[eqcount - 1] = 1;
	S.reverseRun(eqcount - 1);
	return normalize(S.solution);
}

void SLAE::conjugateGradientMethod()
{
	std::vector<double> x;
	x.assign(eqcount, 0);
	std::vector<double> discrep(B - A*x);
	std::vector<double> p(discrep);
	std::vector<double> q;
	double alpha, beta;

	while (norm2(discrep) > eps) {
		q = A * p;
		alpha = scal(discrep, p) / scal(q, p);
		x = x + p * alpha;
		discrep = discrep - q * alpha;
		beta = scal(discrep, q) / scal(p, q);
		p = discrep - p * beta;
	}

	solution = x;

}

double SLAE::maxEigenValue()
{
	double res;

	std::vector<double> y;
	y.assign(eqcount, 1);
	std::vector<double> x(normalize(y));
	std::vector<double> newx;

	std::vector<double> lambda1, lambda2;
	lambda1.assign(eqcount, 0);

	bool notconv(true);
	while (notconv) {
		y = A * x;
		newx = normalize(y);
		lambda2 = y / x;
		notconv = checkConvergence(lambda1, lambda2);
		lambda1 = lambda2;
		x = newx;
	}

	res = average(lambda1);

	return res;
}

double SLAE::minEigenValue()
{
	double res;

	double mEV = maxEigenValue();

	std::vector<std::vector<double> > E(eqcount);
	for (int i(0); i < eqcount; ++i) {
		E[i].assign(varcount, 0);
		E[i][i] = 1;
	}

	std::vector<std::vector<double> > mat(A - mEV * E);

	SLAE S(mat);

	double lambda(S.maxEigenValue());
	res = lambda + mEV;

	return res;
}

