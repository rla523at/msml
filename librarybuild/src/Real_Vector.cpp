#include "Real_Vector.h"
#include "Exception.h"
#include "MSBLAS.h"

#include <iomanip>
#include <cmath>
#include <sstream>


Const_Real_Vector_Wrapper::Const_Real_Vector_Wrapper(const double* coordinate_const_ptr, const int dimension, const int incx)
{
	REQUIRE(0 < dimension, "dimension should greater than 0");
	REQUIRE(coordinate_const_ptr != nullptr, "data ptr should not be nullptr");
	this->dimension_ = dimension;
	this->coordinate_const_ptr_ = coordinate_const_ptr;
	this->inc_ = incx;
}

Const_Real_Vector_Wrapper::Const_Real_Vector_Wrapper(const std::vector<double>& coordinates, const int incx)
{
	const auto dimension = static_cast<int>(std::ceil(coordinates.size() / static_cast<double>(incx)));

	REQUIRE(0 < dimension, "dimension should greater than 0");
	this->dimension_ = dimension;
	this->coordinate_const_ptr_ = coordinates.data();
	this->inc_ = incx;
}


Real_Vector<0> Const_Real_Vector_Wrapper::operator*(const double scalar) const
{
	Real_Vector result(this->dimension_);
	constexpr auto incr = 1;

	ms::blas::cx(result.data(), scalar, this->coordinate_const_ptr_, this->dimension_, incr, this->inc_);

	return result;
}

Real_Vector<0> Const_Real_Vector_Wrapper::operator+(const Const_Real_Vector_Wrapper& other) const
{
	REQUIRE(this->dimension_ == other.dimension_, "other vector should be same size");

	Real_Vector result(this->dimension_);
	constexpr auto incr = 1;

	ms::blas::x_plus_y(result.data(), this->coordinate_const_ptr_, other.coordinate_const_ptr_, this->dimension_, incr, this->inc_, other.inc_);

	return result;
}

Real_Vector<0> Const_Real_Vector_Wrapper::operator-(const Const_Real_Vector_Wrapper& other) const
{
	REQUIRE(this->dimension_ == other.dimension_, "other vector should be same size");

	Real_Vector result(this->dimension_);
	constexpr auto incr = 1;

	ms::blas::x_minus_y(result.data(), this->coordinate_const_ptr_, other.coordinate_const_ptr_, this->dimension_, incr, this->inc_, other.inc_);

	return result;
}

double Const_Real_Vector_Wrapper::operator[](const size_t position) const
{
	REQUIRE(position < this->dimension_, "position should be less then size");
	return this->coordinate_const_ptr_[position * this->inc_];
}

bool Const_Real_Vector_Wrapper::operator==(const Const_Real_Vector_Wrapper& other) const
{
	if (this->dimension_ != other.dimension_)
		return false;

	for (int i = 0; i < this->dimension_; ++i)
	{
		if (this->coordinate_const_ptr_[i*this->inc_] != other.coordinate_const_ptr_[i*other.inc_])
			return false;
	}

	return true;
}

bool Const_Real_Vector_Wrapper::operator!=(const Const_Real_Vector_Wrapper& other) const
{
	return !((*this) == other);
}

double Const_Real_Vector_Wrapper::at(const size_t position) const
{
	REQUIRE(position < this->dimension_, "position can't exceed given range");
	return this->coordinate_const_ptr_[position * this->inc_];
}

const double* Const_Real_Vector_Wrapper::begin(void) const
{
	return this->coordinate_const_ptr_;
}

const double* Const_Real_Vector_Wrapper::data(void) const
{
	return this->coordinate_const_ptr_;
}

const double* Const_Real_Vector_Wrapper::end(void) const
{
	const auto end_index = this->inc_ * (this->dimension_ - 1) + 1;
	return this->coordinate_const_ptr_ + end_index;
}

double Const_Real_Vector_Wrapper::L1_norm(void) const
{
	return ms::blas::abs_sum_x(this->coordinate_const_ptr_, this->dimension_, this->inc_);
}

double Const_Real_Vector_Wrapper::L2_norm(void) const
{
	return std::sqrt(this->inner_product(*this));
}

double Const_Real_Vector_Wrapper::Linf_norm(void) const
{
	const auto pos = ms::blas::find_maximum_element_pos(this->coordinate_const_ptr_, this->dimension_, this->inc_);
	return this->at(pos);
}

double Const_Real_Vector_Wrapper::inner_product(const Const_Real_Vector_Wrapper& other) const
{
	REQUIRE(this->dimension_ == other.dimension_, "two vector should be same size");

	return ms::blas::x_dot_y(this->coordinate_const_ptr_, other.coordinate_const_ptr_, this->dimension_);
}

size_t Const_Real_Vector_Wrapper::dimension(void) const
{
	return this->dimension_;
}

Const_Real_Vector_Wrapper Const_Real_Vector_Wrapper::part(const int start_pos, const int end_pos) const
{
	REQUIRE(start_pos < end_pos, "start position should be smaller than end position");
	REQUIRE(end_pos <= this->dimension_, "end position should be smaller or equal to dimension");

	const auto start_index = start_pos * this->inc_;
	const auto start_ptr = this->coordinate_const_ptr_ + start_index;
	const auto dimension_of_part = end_pos - start_pos;

	return { start_ptr, dimension_of_part, this->inc_ };
}

std::string Const_Real_Vector_Wrapper::to_string(void) const
{
	std::ostringstream oss;
	oss << std::setprecision(16) << std::showpoint << std::left;
	for (int i = 0; i < this->dimension_; ++i)
	{
		oss << std::setw(25) << this->coordinate_const_ptr_[i];
	}
	return oss.str();
}

template<>
Real_Vector<3> Const_Real_Vector_Wrapper::cross_product<2>(const Const_Real_Vector_Wrapper& other) const
{
	Real_Vector<3> result;
	result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

	return result;
}

template<>
Real_Vector<3> Const_Real_Vector_Wrapper::cross_product<3>(const Const_Real_Vector_Wrapper& other) const
{
	Real_Vector<3> result;
	result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
	result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
	result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

	return result;
}

void Real_Vector_Wrapper::operator*=(const double constant)
{
	ms::blas::cx(constant, this->coordinate_ptr_, this->dimension_, this->inc_);
}

void Real_Vector_Wrapper::operator+=(const Const_Real_Vector_Wrapper& other)
{
	REQUIRE(this->dimension_ == other.dimension(), "other vector should be same size");
	ms::blas::x_plus_assign_y(this->coordinate_ptr_, other.data(), this->dimension_, this->inc_, other.inc());
}

void Real_Vector_Wrapper::operator-=(const Const_Real_Vector_Wrapper& other)
{
	REQUIRE(this->dimension_ == other.dimension(), "other vector should be same size");
	ms::blas::x_minus_assign_y(this->coordinate_ptr_, other.data(), this->dimension_, this->inc_, other.inc());
}

double& Real_Vector_Wrapper::operator[](const size_t position)
{
	REQUIRE(position < this->dimension_, "position should be less then size");
	return this->coordinate_ptr_[position];
}

double& Real_Vector_Wrapper::at(const size_t position)
{
	REQUIRE(position < this->dimension_, "position should be less then size");
	return this->coordinate_ptr_[position];
}

void Real_Vector_Wrapper::change_value(const Const_Real_Vector_Wrapper& other)
{
	REQUIRE(this->dimension_ == other.dimension(), "other vector should be same size");
	ms::blas::copy(this->coordinate_ptr_, other.data(), this->dimension_, this->inc_, other.inc());
}


double* Real_Vector_Wrapper::data(void)
{
	return this->coordinate_ptr_;
}

void Real_Vector_Wrapper::initalize(void)
{
	for (int i = 0; i < this->dimension_; ++i)
	{
		this->coordinate_ptr_[i] = 0.0;
	}
}

void Real_Vector_Wrapper::normalize(void)
{
	*this *= 1.0 / this->L2_norm();
}

Real_Vector<0>::Real_Vector(const int dimension)
	: coordinates_(dimension)
{
	this->dimension_ = static_cast<int>(this->coordinates_.size());
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = 1;
}

Real_Vector<0>::Real_Vector(const std::initializer_list<double> list)
	:coordinates_(list)
{
	this->dimension_ = static_cast<int>(this->coordinates_.size());
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = 1;
}

Real_Vector<0>::Real_Vector(const std::vector<double>& values)
	:coordinates_(values)
{
	this->dimension_ = static_cast<int>(this->coordinates_.size());
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = 1;
}

Real_Vector<0>::Real_Vector(std::vector<double>&& values)
	:coordinates_(std::move(values))
{
	this->dimension_ = static_cast<int>(this->coordinates_.size());
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = 1;
};

Real_Vector<0>::Real_Vector(const Real_Vector& other)
	:coordinates_(other.coordinates_)
{
	this->dimension_ = other.dimension_;
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = other.inc_;
}

Real_Vector<0>::Real_Vector(Real_Vector&& other) noexcept
	:coordinates_(std::move(other.coordinates_))
{
	this->dimension_ = other.dimension_;
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = other.inc_;

	other.dimension_ = 0;
	other.coordinate_const_ptr_ = nullptr;
	other.coordinate_ptr_ = nullptr;
}

void Real_Vector<0>::operator=(const Real_Vector& other)
{
	this->coordinates_ = other.coordinates_;
	this->dimension_ = other.dimension_;
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = other.inc_;
}

void Real_Vector<0>::operator=(Real_Vector&& other) noexcept
{
	this->coordinates_ = std::move(other.coordinates_);
	this->dimension_ = other.dimension_;
	this->coordinate_const_ptr_ = this->coordinates_.data();
	this->coordinate_ptr_ = this->coordinates_.data();
	this->inc_ = other.inc_;
}

Real_Vector<0> operator*(const double constant, const Const_Real_Vector_Wrapper& x)
{
	return x * constant;
}