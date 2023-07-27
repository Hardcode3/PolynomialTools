#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <eigen/Eigen/Dense>

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

class Polynomial
{
public:
	Polynomial() = delete;
	Polynomial(std::vector<double> coefficients);
	Polynomial(const Polynomial& polynomial);
	~Polynomial() { ; }

	double evaluate(const double& x);
	void normalize();
	Eigen::VectorXcd compute_roots();

	Eigen::MatrixXd get_companion_matrix() const;
	const std::vector<double>& get_coefficients() const { return coefficients_; }
	const double& get_coefficient(const std::size_t& index) const;
	const Eigen::VectorXcd& get_roots() const { return roots_; }
	std::size_t get_first_non_null_coeff() const { return first_non_null_coef_index_; }
	std::size_t order() const { return order_; }
	std::size_t get_null_coef_nb() const { return null_coeff_nb_; }

	void print() const { std::cout << *this << std::endl; }
	friend std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial);
	const double& operator[](std::size_t index) const;
	double& operator[](std::size_t index);
	Polynomial& operator=(const Polynomial& other);
	bool is_approx(const Polynomial& other, const double tolerance = 1e-10);
	bool operator==(const Polynomial& other);
	bool operator!=(const Polynomial& other);

	const double TOLERANCE = 1e-10;
	const int ND_VALUE = -99999;

private:
	std::vector<double> coefficients_;
	Eigen::VectorXcd roots_;
	std::size_t first_non_null_coef_index_ = 0;
	std::size_t order_ = 0;
	std::size_t null_coeff_nb_ = 0;
};

#endif