#include "polynomial.hpp"

Polynomial::Polynomial(std::vector<double> coefficients)
{
	if (coefficients.size() == 0)
		throw std::invalid_argument("There is no coefficient for this polynomial");
	coefficients_ = coefficients;

	order_ = coefficients_.size() - 1;

	// check if the polynomial is null
	for (double& coef : coefficients_)
	{
		if (coef == 0)
			null_coeff_nb_++;
	}
	if (null_coeff_nb_ == order_ + 1)
		throw std::invalid_argument("Can not instanciate a null polynomial (no operation to do on it)");

	// find the first non null coefficient index
	while (coefficients_[first_non_null_coef_index_] == 0)
	{
		first_non_null_coef_index_++;
	}

}

Polynomial::Polynomial(const Polynomial& other)
{
	this->coefficients_ = other.coefficients_;
	this->roots_ = other.roots_;
	this->first_non_null_coef_index_ = other.first_non_null_coef_index_;
	this->order_ = other.order_;
	this->null_coeff_nb_ = other.null_coeff_nb_;
}

double Polynomial::evaluate(const double& x)
{
	double value = 0;
	for (int i = 0; i <= order_; i++)
		value += coefficients_[order_ - i] * std::pow(x, i);
	return value;
}

void Polynomial::normalize()
{
	double first_non_null_coeff = coefficients_[first_non_null_coef_index_];
	for (double& element : coefficients_)
		element /= first_non_null_coeff;
}

Eigen::VectorXcd Polynomial::compute_roots()
{
	if (order_ == 0)
		return Eigen::VectorXcd();

	// construct the companion matrix for the polynomial (dynamic size)
	Eigen::MatrixXd companion_matrix(order_, order_);
	companion_matrix.setZero();

	for (std::size_t i = 0; i < order_ - 1; i++)
		companion_matrix(i + 1, i) = 1;

	for (std::size_t i = 0; i < order_; i++)
		companion_matrix(0, i) = -coefficients_[i+1] / coefficients_[first_non_null_coef_index_];

	// solve for the eigen values that correspond to the roots of the monic polynomial
	Eigen::EigenSolver<Eigen::MatrixXd> solver(companion_matrix);
	Eigen::VectorXcd roots_= solver.eigenvalues();

	return roots_;
}

Eigen::MatrixXd Polynomial::get_companion_matrix() const
{
	Eigen::MatrixXd companion_matrix(order_, order_);
	companion_matrix.setZero();

	for (std::size_t i = 0; i < order_ - 1; i++)
		companion_matrix(i + 1, i) = 1;

	for (std::size_t i = 0; i < order_; i++)
		companion_matrix(0, i) = -coefficients_[i + 1] / coefficients_[0];

	return companion_matrix;
}

const double& Polynomial::get_coefficient(const std::size_t& index) const
{
	if (0 <= index && index <= order_)
		return coefficients_[index];
	else
		throw std::invalid_argument("The specified index is out of bounds");
}

std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial)
{
	os << std::setprecision(15);
	os << "| Polynomial object at adress " << &polynomial << "\n";
	os << "| Polynomial order: " << polynomial.order_ << "\n";
	os << "| First non null coefficient index: " << polynomial.first_non_null_coef_index_ << "\n";
	os << "| Number of null coefficients:  " << polynomial.null_coeff_nb_ << "\n";
	os << "| P(x) = ";
	for (std::size_t i = 0; i < polynomial.order_ ; i++)
	{
		os << polynomial.coefficients_[i] << " * x^" << polynomial.order_ - i << "  +  ";
	}
	os << polynomial.coefficients_[polynomial.order_] << "\n";
	return os;
}

double& Polynomial::operator[](std::size_t index)
{
	if (0 <= index && index <= order_) \
	{
		return coefficients_[index];
	}
	else
	{
		throw std::out_of_range("Index out of range");
	}
}

const double& Polynomial::operator[](std::size_t index) const {
	if (0 <= index && index <= order_) \
	{
		return coefficients_[index];
	}
	else
	{
		throw std::out_of_range("Index out of range");
	}
}

Polynomial& Polynomial::operator=(const Polynomial& other) 
{
	if (this == &other)
		return *this;

	this->coefficients_ = other.coefficients_;
	this->roots_ = other.roots_;
	this->first_non_null_coef_index_ = other.first_non_null_coef_index_;
	this->order_ = other.first_non_null_coef_index_;
	this->null_coeff_nb_ = other.null_coeff_nb_;

	// Return a reference to the current object
	return *this;
}

bool Polynomial::is_approx(const Polynomial& other, const double tolerance)
{
	if (this == &other)
		return true;

	if (this->order_ != other.order_)
		return false;

	for (std::size_t i = 0; i < this->order_; i++)
		if (this->coefficients_[i] - other.coefficients_[i] > tolerance)
			return false;

	return true;
}

bool Polynomial::operator==(const Polynomial& other)
{
	if (this == &other)
		return true;

	if (this->order_ != other.order_)
		return false;

	for (std::size_t i = 0; i < this->order_; i++)
		if (this->coefficients_[i] - other.coefficients_[i] > TOLERANCE)
			return false;
		
	return true;
}

bool Polynomial::operator!=(const Polynomial& other)
{
	return !(*this==other);
}


