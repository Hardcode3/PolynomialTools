#include "polynomial.hpp"

int test_normalisation(const double& tol = 1e-7)
{
	// 1e7 tolerance to compare the double precision numbers with the hardcoded expected values
	Polynomial p({ 6, 2, 3, 5 });
	Polynomial p_monic({ 1, 0.333333333333333, 0.5, 0.8333333333333334 });
	
	if (p == p_monic)
		return 1;
	if (p.is_approx(p_monic, tol))
		return 1;

	p.normalize();

	if (p != p_monic)
		return 1;

	// check normalisation
	for (int i = 0; i < p.order(); i++)
		if (p[i] - p_monic[i] > p.TOLERANCE)
			return 1;
	
	Polynomial p_cpy(p);

	if (p != p_cpy)
		return 1;

	Eigen::VectorXcd roots = p.compute_roots();
	Eigen::VectorXcd roots_expected(p.order());
	roots_expected <<	std::complex<double>(0.26653043, 0.94382327),
						std::complex<double>(0.26653043, -0.94382327),
						std::complex<double>(-0.8663942, 0.);
	
	if (!roots.isApprox(roots_expected, tol))
		return 1;

	Polynomial p_1({ 6.9, 1.9, 0.03, 0.15 });
	roots = p_1.compute_roots();
	roots_expected <<	std::complex<double>(-0.40021937, 0.),
						std::complex<double>(0.06242853, 0.22454558),
						std::complex<double>(0.06242853, -0.22454558);
	if (p_1.evaluate(0) - 0.15 > tol)
		return 1;

	if (p_1.evaluate(1.2) - 14.8452> tol)
		return 1;

	if (!roots.isApprox(roots_expected, tol))
		return 1;

	Polynomial p_2({ 6.9, 1.9 });
	roots = p_2.compute_roots();
	roots_expected = Eigen::VectorXcd(p_2.order());
	roots_expected << std::complex<double>(-0.27536232, 0.);

	if (!roots.isApprox(roots_expected, tol))
		return 1;

	Polynomial p_3({ 6.9 });
	std::cout << p_3.compute_roots();
	roots = p_3.compute_roots();

	if (roots.size() != 0)
		return 1;
	
	Polynomial p_4({ 1.5, -1.9, 0.03, 0.15, -4.534, 0.0023, 2.9842, -1.545 });
	roots = Eigen::VectorXcd(p_4.compute_roots());
	roots_expected = Eigen::VectorXcd(p_4.order());
	roots_expected <<	std::complex<double>(1.72155401, 0.),
						std::complex<double>(0.22626978, 1.3097414),
						std::complex<double>(0.22626978, -1.3097414),
						std::complex<double>(-0.95227588, 0.32329302),
						std::complex<double>(-0.95227588, -0.32329302),
						std::complex<double>(0.49856243, +0.29377346),
						std::complex<double>(0.49856243, -0.29377346);

	if (!roots.isApprox(roots_expected, tol))
		return 1;

	return 0;
}

int main()
{
	if (test_normalisation() == 1)
		return 1;
	return 0;
}