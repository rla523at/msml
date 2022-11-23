#pragma once
#include <array>
#include <string>
#include <vector>

#define COMPILE_TIME_REQUIREMENT static_assert

namespace ms
{
	template <typename... Args>
	inline constexpr bool are_arithmetics = (... && std::is_arithmetic_v<Args>);

	template <typename T>
	using if_is_iterator = std::enable_if_t<std::_Is_iterator_v<T>, bool>;

	template <typename... Args>
	using if_is_more_than_two = std::enable_if_t<2 <= sizeof...(Args), bool>;
}

template <int dim>
class Real_Vector;

class Const_Real_Vector_Wrapper
{
public:
	Const_Real_Vector_Wrapper(void) = default;
	Const_Real_Vector_Wrapper(const double* coordinate_ptr, const int dimension, const int incx = 1);
	Const_Real_Vector_Wrapper(const std::vector<double>& coordinates, const int incx = 1);

public://Query
	Real_Vector<0> operator*(const double scalar) const;
	Real_Vector<0> operator-(const Const_Real_Vector_Wrapper& other) const;
	Real_Vector<0> operator+(const Const_Real_Vector_Wrapper& other) const;
	double operator[](const size_t position) const;
	bool operator==(const Const_Real_Vector_Wrapper& other) const;
	bool operator!=(const Const_Real_Vector_Wrapper& other) const;

	double at(const size_t position) const;
	const double* begin(void) const;
	const double* data(void) const;
	size_t dimension(void) const;
	const double* end(void) const;
	double inner_product(const Const_Real_Vector_Wrapper& other) const; 
	int inc(void) const { return this->inc_; };
	double L1_norm(void) const;
	double L2_norm(void) const;
	double Linf_norm(void) const;	
	// [start_pos, end_pos)
	Const_Real_Vector_Wrapper part(const int start_pos, const int end_pos) const; 
	std::string to_string(void) const;

	template<int dim>
	Real_Vector<3> cross_product(const Const_Real_Vector_Wrapper& other) const;

protected:
	int dimension_ = 0;
	const double* coordinate_const_ptr_ = nullptr;
	int inc_ = 0;
};

class Real_Vector_Wrapper : public Const_Real_Vector_Wrapper
{
public:
	Real_Vector_Wrapper(void) = default;
	Real_Vector_Wrapper(double* coordinate_ptr, const int dimension, const int incx = 1)
		: Const_Real_Vector_Wrapper(coordinate_ptr, dimension, incx)
		, coordinate_ptr_(coordinate_ptr) {};
	Real_Vector_Wrapper(std::vector<double>& coordinates, const int incx = 1)
		: Const_Real_Vector_Wrapper(coordinates, incx)
		, coordinate_ptr_(coordinates.data()) {};

public://Command
	void operator*=(const double constant);
	void operator+=(const Const_Real_Vector_Wrapper& other);
	void operator-=(const Const_Real_Vector_Wrapper& other);
	double& operator[](const size_t position);

	double& at(const size_t position);	
	void change_value(const Const_Real_Vector_Wrapper& other);
	double* data(void);
	void initalize(void);
	void normalize(void);

	//resolve base class method hiding problem (overloading across scope problem)
	using Const_Real_Vector_Wrapper::operator[];
	using Const_Real_Vector_Wrapper::at;
	using Const_Real_Vector_Wrapper::data;

protected:
	double* coordinate_ptr_ = nullptr;
};

template <int Dim = 0>
class Real_Vector : public Real_Vector_Wrapper
{
	COMPILE_TIME_REQUIREMENT(0 <= Dim, "dimension shold be positive");

public:
	Real_Vector(void)
		: Real_Vector_Wrapper(this->coordinates_.data(), Dim) {};
	Real_Vector(const Real_Vector& other)
		: Real_Vector_Wrapper(this->coordinates_.data(), Dim)
		, coordinates_{ other.coordinates_ } {};
	template <typename ... Args>
	Real_Vector(const Args&... args)
		: Real_Vector_Wrapper(this->coordinates_.data(), Dim)
		, coordinates_{ static_cast<double>(args)... }
	{
		COMPILE_TIME_REQUIREMENT(sizeof...(Args) <= Dim, "Number of arguments can not exceed dimension");
		COMPILE_TIME_REQUIREMENT(ms::are_arithmetics<Args...>, "Every arguments should be arithmetics");
	}

public:
	void operator=(const Real_Vector& other)
	{
		this->coordinates_ = other.coordinates_;
	}

private:
	std::array<double, Dim> coordinates_ = { 0 };
};

//user-defined deduction guides
template <typename Iter, ms::if_is_iterator<Iter> = true> 
Real_Vector(Iter iter1, Iter iter2)->Real_Vector<0>;

template <typename... Args, ms::if_is_more_than_two<Args...> = true> 
Real_Vector(Args... args)->Real_Vector<sizeof...(Args)>;


template <>
class Real_Vector<0> : public Real_Vector_Wrapper
{
public:
	Real_Vector(void) = default;
	explicit Real_Vector(const int dimension);
	Real_Vector(const std::initializer_list<double> list);
	Real_Vector(const std::vector<double>&values);
	Real_Vector(std::vector<double>&& values);
	Real_Vector(const Real_Vector& other);
	Real_Vector(Real_Vector && other) noexcept;

	template <typename Iter, ms::if_is_iterator<Iter> = true>
	Real_Vector(Iter first, Iter last)
		: coordinates_(first, last)
	{
		this->dimension_ = static_cast<int>(this->coordinates_.size());
		this->coordinate_const_ptr_ = this->coordinates_.data();
		this->coordinate_ptr_ = this->coordinates_.data();
		this->inc_ = 1;
	};

public://Command	
	void operator=(const Real_Vector & other);
	void operator=(Real_Vector && other) noexcept;

private:
	std::vector<double> coordinates_;
};

Real_Vector<0> operator*(const double constant, const Const_Real_Vector_Wrapper& x);