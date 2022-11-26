#pragma once
#include "test_Matrix_API.h"
#include "../MS_Math_Library/src/Real_Vector.h"

//for google test cout message 
std::ostream& operator<<(std::ostream& os, const Matrix_Wrapper& mw)
{
	return os << mw.to_string();
}
//for google test cout message 


TEST(Matrix_Wrapper, at_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	std::vector<double> val = { 1,2,3,4,5,6,7,8,9 };
	constexpr auto leading_dimension = 3;
	
	Matrix_Wrapper m = { num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1 2 ��T	 =	��	1 4 7	��
	//			�� 4 5 ��				��	2 5 8	��
	//			�� 7 8 ��

	m.at(1, 1) = 0;

	std::vector<double> ref_val = { 1,2,3,4,0,6,7,8,9 };
	EXPECT_EQ(val, ref_val);
}
TEST(Matrix_Wrapper, operator_addition_1)
{
	constexpr auto num_rows1 = 2;
	constexpr auto num_columns1 = 3;
	constexpr auto leading_dimension1 = 5;
	std::vector<double> val1 = { 1,2,3,4,5,6,7,8,9,10 };
	Matrix_Wrapper mw = { num_rows1, num_columns1, val1.data(), leading_dimension1 };
	// cmw1 =	�� 1	2	3 ��
	//				�� 6	7	8 ��

	constexpr auto num_rows2 = 3;
	constexpr auto num_columns2 = 2;
	constexpr auto leading_dimension2 = 3;
	std::vector<double> val2 = { 1,2,3,4,5,6,7,8,9 };
	const Matrix_Wrapper cmw = { num_rows2, num_columns2, val2.data(), leading_dimension2, Transpose_Type::transposed };
	// cmw2 =	�� 1	2	��T	=	�� 1	4	7 ��		
	//				�� 4	5	��			�� 2	5	8 ��
	//				�� 7	8	��

	mw += cmw;
	std::vector<double> ref = { 2,6,10,4,5,8,12,16,9,10 };
	EXPECT_EQ(val1, ref);
}
TEST(Matrix_Wrapper, operator_assign_scalar_multiplication_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	std::vector<double> val = { 1,2,3,4,5,6,7,8,9 };
	constexpr auto leading_dimension = 3;
	
	Matrix_Wrapper m = { num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1 2 ��T	 =	��	1 4 7	��
	//			�� 4 5 ��				��	2 5 8	��
	//			�� 7 8 ��

	m *= 3;

	std::vector<double> ref_val = { 3,6,3,12,15,6,21,24,9 };
	EXPECT_EQ(val, ref_val);
}
TEST(Matrix_Wrapper, operator_assign_scalar_multiplication_2)
{
	constexpr int num_rows = 10000;
	constexpr int num_columns = 10000;
	constexpr int leading_dimension = num_columns;
	constexpr double initial_value = 1.0101;

	std::vector<double> values(num_rows * leading_dimension, initial_value);

	Matrix_Wrapper m = { num_rows, num_columns, values.data(), leading_dimension};

	constexpr double c = -1.0;
	m *= c;

	std::vector<double> ref_values( num_rows * leading_dimension, -initial_value);
	EXPECT_EQ(values, ref_values);
}
TEST(Matrix_Wrapper, column_vector_wrapper_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	std::vector<double> val = { 1,2,3,4,5,6,7,8,9 };
	constexpr auto leading_dimension = 3;
	
	Matrix_Wrapper m = { num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1 2 ��T	 =	��	1 4 7	��
	//			�� 4 5 ��				��	2 5 8	��
	//			�� 7 8 ��

	constexpr auto column_index = 0;
	auto column_vector = m.column_vector_wrapper(column_index);

	std::vector<double> new_val = { 5,6 };
	column_vector.change_value(new_val);

	std::vector<double> ref_val = { 5,6,3,4,5,6,7,8,9 };
	EXPECT_EQ(val, ref_val);
}
TEST(Matrix_Wrapper, column_vector_wrapper_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	std::vector<double> val = { 1,2,3,4,5,6,7,8,9 };
	constexpr auto leading_dimension = 3;
	
	Matrix_Wrapper m = { num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1 2 ��T	 =	��	1 4 7	��
	//			�� 4 5 ��				��	2 5 8	��
	//			�� 7 8 ��

	constexpr auto column_index = 2;
	auto column_vector = m.column_vector_wrapper(column_index);

	std::vector<double> new_val = { 5,6 };
	column_vector.change_value(new_val);

	std::vector<double> ref_val = { 1,2,3,4,5,6,5,6,9 };
	EXPECT_EQ(val, ref_val);
}
TEST(Matrix_Wrapper, change_values_1)
{
	constexpr auto num_rows1 = 3;
	constexpr auto num_columns1 = 2;
	std::vector<double> vals1 = { 1,2,3,4,5,6,7,8,9 };
	constexpr auto leading_dimension1 = 3;

	Matrix_Wrapper m1 = { num_rows1, num_columns1, vals1.data(), leading_dimension1, Transpose_Type::transposed };
	// m1 =	�� 1 2 ��T	 =	��	1 4 7	��
	//			�� 4 5 ��				��	2 5 8	��
	//			�� 7 8 ��

	constexpr auto num_rows2 = 2;
	constexpr auto num_columns2 = 3;
	std::vector<double> vals2 = { 1,2,3,4,5,6,7,8,9,10 };
	constexpr auto leading_dimension2 = 5;

	Const_Matrix_Wrapper m2 = { num_rows2, num_columns2, vals2.data(), leading_dimension2, Transpose_Type::not_transposed };
	// m2 =	�� 1  2  3 ��
	//			�� 6  7  8 ��

	m1.change_value(m2);

	std::vector<double> ref_vals = { 1,6,3,2,7,6,3,8,9 };
	EXPECT_EQ(vals1, ref_vals);
}
TEST(Matrix_Wrapper, change_values_2)
{
	constexpr auto num_rows1 = 4;
	constexpr auto num_columns1 = 3;
	std::vector<double> vals1 = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 };
	constexpr auto leading_dimension1 = 4;

	Matrix_Wrapper m1 = { num_rows1, num_columns1, vals1.data(), leading_dimension1, Transpose_Type::not_transposed };
	// m1 =	�� 1  2  3  ��
	//			�� 5  6  7  ��
	//			�� 9  10 11 ��
	//			�� 13 14 15 ��

	constexpr auto num_rows2 = 3;
	constexpr auto num_columns2 = 4;
	std::vector<double> vals2 = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
	constexpr auto leading_dimension = 5;

	Const_Matrix_Wrapper m2 = { num_rows2, num_columns2, vals2.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1  2  3  4  ��T	 =	��	1 6 11 ��
	//			�� 6  7  8  9  ��				��	2 7 12 ��
	//			�� 11 12 13 14	��				�� 3	8 13 ��
	//														�� 4 9 14 ��

	m1.change_value(m2);

	std::vector<double> ref_vals = { 1,6,11,4,2,7,12,8,3,8,13,12,4,9,14,16 };
	EXPECT_EQ(vals1, ref_vals);
}
TEST(Matrix_Wrapper, row_vector_wrapper_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	std::vector<double> val = { 1,2,3,4,5,6,7,8,9 };
	constexpr auto leading_dimension = 3;
	
	Matrix_Wrapper m = { num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1 2 ��T	 =	��	1 4 7	��
	//			�� 4 5 ��				��	2 5 8	��
	//			�� 7 8 ��

	constexpr auto row_index = 0;
	auto row_vector = m.row_vector_wrapper(row_index);

	std::vector<double> new_val = { 1,2,3 };
	row_vector.change_value(new_val);

	std::vector<double> ref_val = { 1,2,3,2,5,6,3,8,9 };
	EXPECT_EQ(val, ref_val);
}
TEST(Matrix_Wrapper, part_1)
{
	constexpr auto num_row = 2;
	constexpr auto num_column = 4;
	std::vector<double> val = { 1,2,3,4,5,6,7,8 };
	Matrix_Wrapper m = { num_row, num_column, val.data() };
	// m =	�� 1 2 3 4 ��	
	//			�� 5 6 7 8 ��		

	constexpr auto start_row_index = 0;
	constexpr auto end_row_index = 2;
	constexpr auto start_column_index = 1;
	constexpr auto end_column_index = 3;
	auto result = m.part(start_row_index, end_row_index, start_column_index, end_column_index);

	result *= -1;
	std::vector<double> ref_val = { 1, -2, -3, 4, 5 ,-6, -7, 8 };

	EXPECT_EQ(val, ref_val);
}
TEST(Matrix_Wrapper, part_2)
{
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 4;
	std::vector<double> vals = { 1,2,3,4,5,6,7,8 };
	Matrix_Wrapper m = { num_rows, num_columns, vals.data() };
	// m =	�� 1 2 3 4 ��	
	//			�� 5 6 7 8 ��		

	constexpr auto start_row_index = 0;
	constexpr auto end_row_index = 2;
	constexpr auto start_column_index = 0;
	constexpr auto end_column_index = 3;
	auto result = m.part(start_row_index, end_row_index, start_column_index, end_column_index);

	result *= -1.0;

	std::vector<double> ref_vals = { -1, -2, -3, 4, -5 ,-6, -7, 8 };
	EXPECT_EQ(vals, ref_vals);
}
TEST(Matrix_Wrapper, part_3)
{
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 4;
	std::vector<double> vals = { 1,2,3,4,5,6,7,8 };
	
	Matrix_Wrapper m = { num_rows, num_columns, vals.data(), Transpose_Type::transposed };
	// m =	�� 1 2 3 4 ��T	 =	��	1 5 ��
	//			�� 5 6 7 8 ��				��	2 6 ��
	//												�� 3 7 ��
	//												�� 4 8 ��

	constexpr auto start_row_index = 1;
	constexpr auto end_row_index = 3;
	constexpr auto start_column_index = 0;
	constexpr auto end_column_index = 2;
	auto result = m.part(start_row_index, end_row_index, start_column_index, end_column_index);

	result *= -1.0;

	std::vector<double> ref_vals = { 1,-2,-3,4,5,-6,-7,8 };
	EXPECT_EQ(vals, ref_vals);
}
TEST(Matrix_Wrapper, part_4)
{
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 4;
	std::vector<double> vals = { 1,2,3,4,5,6,7,8 };
	
	Matrix_Wrapper m = { num_rows, num_columns, vals.data(), Transpose_Type::transposed };
	// m =	�� 1 2 3 4 ��T	 =	��	1 5 ��
	//			�� 5 6 7 8 ��				��	2 6 ��
	//												�� 3 7 ��
	//												�� 4 8 ��

	constexpr auto start_row_index = 1;
	constexpr auto end_row_index = 4;
	constexpr auto start_column_index = 0;
	constexpr auto end_column_index = 2;
	auto result = m.part(start_row_index, end_row_index, start_column_index, end_column_index);

	result *= -1.0;

	std::vector<double> ref_vals = { 1,-2,-3,-4,5,-6,-7,-8 };
	EXPECT_EQ(vals, ref_vals);
}
TEST(Matrix_Wrapper, part_5)
{
	constexpr auto num_rows = 3;
	constexpr auto num_columns = 4;
	std::vector<double> vals = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
	constexpr auto leading_dimension = 5;
	
	Matrix_Wrapper m = { num_rows, num_columns, vals.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1  2  3  4  ��T	 =	��	1 6 11 ��
	//			�� 6  7  8  9  ��				��	2 7 12 ��
	//			�� 11 12 13 14	��				�� 3	8 13 ��
	//														�� 4 9 14 ��

	constexpr auto start_row_index = 0;
	constexpr auto end_row_index = 3;
	constexpr auto start_column_index = 1;
	constexpr auto end_column_index = 3;
	auto result = m.part(start_row_index, end_row_index, start_column_index, end_column_index);

	result *= -1.0;

	std::vector<double> ref_vals = { 1,2,3,4,5,-6,-7,-8,9,10,-11,-12,-13,14,15 };
	EXPECT_EQ(vals, ref_vals);
}
TEST(Matrix_Wrapper, part_6)
{
	constexpr auto num_rows = 3;
	constexpr auto num_columns = 4;
	std::vector<double> vals = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
	constexpr auto leading_dimension = 5;

	Matrix_Wrapper m = { num_rows, num_columns, vals.data(), leading_dimension, Transpose_Type::transposed };
	// m =	�� 1  2  3  4  ��T	 =	��	1 6 11 ��
	//			�� 6  7  8  9  ��				��	2 7 12 ��
	//			�� 11 12 13 14	��				�� 3	8 13 ��
	//														�� 4 9 14 ��

	constexpr auto start_row_index = 2;
	constexpr auto end_row_index = 4;
	constexpr auto start_column_index = 0;
	constexpr auto end_column_index = 3;
	auto result = m.part(start_row_index, end_row_index, start_column_index, end_column_index);

	result *= -1.0;

	std::vector<double> ref_vals = { 1,2,-3,-4,5,6,7,-8,-9,10,11,12,-13,-14,15 };
	EXPECT_EQ(vals, ref_vals);
}

TEST(Matrix_Wrapper, inverse_1)
{
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 2;
	std::vector<double> val = { 2,1,1,1 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data() };
	mw.inverse();

	std::vector<double> ref_val = { 1,-1,-1,2 };
	Matrix_Wrapper ref = { num_rows, num_columns, ref_val.data() };
	EXPECT_EQ(mw, ref);
}
TEST(Matrix_Wrapper, inverse_2)
{
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 2;
	std::vector<double> val = { 1,2,3,4 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data() };
	mw.inverse();

	std::vector<double> ref_val = { -2,1,1.5,-0.5 };
	Matrix_Wrapper ref = { num_rows, num_columns, ref_val.data() };
	
	//EXPECT_EQ(mw, ref); machine error occur!
	matrix_test_API::compare_considering_4ULP(mw, ref); // unit in the last place
} 
TEST(Matrix_Wrapper, inverse_3)
{
	constexpr auto num_rows = 5;
	constexpr auto num_columns = 5;
	std::vector<double> val = { 1, 2, 3, 4, 5, 2, 13, 4, 11, 6, 1, 6, 5, 4, 3, 4, 1, 12, 7, 8, 9, 8, 7, 6, 5 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data() };
	mw.inverse();

	std::vector<double> ref_val = { -0.03333333333333333, 0, -0.1666666666666667, 0,   0.1333333333333333, 0.09375, -0.0625, 0.25, -0.125, 0.03125, -0.1541666666666667, -0.04166666666666667, 0.1666666666666667, 0.08333333333333333, -0.02916666666666667, -0.3395833333333333, 0.2708333333333333, -0.4166666666666667, 0.2083333333333333, -0.06875, 0.5333333333333333, -0.1666666666666667, 0.1666666666666667, -0.1666666666666667, 0.03333333333333333 };
	Matrix_Wrapper ref = { num_rows, num_columns, ref_val.data() };

	//error exceed 4ULP
	const auto epsilon = 1.0E-15;
	matrix_test_API::compare_considering_epsilon(mw, ref, epsilon);
}
TEST(Matrix_Wrapper, inverse_4)
{
	constexpr auto num_rows = 5;
	constexpr auto num_columns = 5;
	std::vector<double> val = { 1, 2, 3, 4, 5,  2, 3, 4, 1, 6, 1, 6, 5, 4, 3, 4, 1, 2, 7, 8, 9, 8, 7, 6, 5 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data() };
	mw.inverse();

	std::vector<double> ref_val = { -0.033333333333333,  -0.000000000000000,  -0.166666666666667,  -0.000000000000000,   0.133333333333333,  -1.291666666666667,   0.250000000000000,   0.666666666666667,   0.500000000000000,  -0.208333333333333,   1.616666666666666,  -0.250000000000000,  -0.666666666666667,  -0.750000000000000,   0.283333333333333,   0.275000000000000,  -0.250000000000000,   0.000000000000000,   0.000000000000000,   0.025000000000000,  -0.466666666666667,   0.250000000000000,   0.166666666666667,   0.250000000000000,  -0.133333333333333 };
	Matrix_Wrapper ref = { num_rows, num_columns, ref_val.data() };

	const auto epsilon = 1.0E-15;
	matrix_test_API::compare_considering_epsilon(mw, ref, epsilon);
}
TEST(Matrix_Wrapper, inverse_5)
{
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> val = { 2,1,3,1,1,3 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data(), leading_dimension };
	mw.inverse();

	std::vector<double> ref_val = { 1,-1,3, -1,2,3 };
	EXPECT_EQ(val, ref_val);
}

TEST(Matrix, transpose_1)
{
	constexpr auto num_rows = 3;
	constexpr auto num_columns = 2;
	std::vector<double> val = { 1,2,3,4,5,6 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data() };
	mw.transpose();

	constexpr auto ref_num_rows = 2;
	constexpr auto ref_num_columns = 3;
	std::vector<double> ref_val = { 1,3,5,2,4,6 };
	Matrix_Wrapper ref = { ref_num_rows, ref_num_columns, ref_val.data() };

	EXPECT_EQ(mw, ref);
}
TEST(Matrix, transpose_2)
{
	constexpr auto num_rows = 3;
	constexpr auto num_columns = 2;
	constexpr auto leading_diemnsion = 3;
	std::vector<double> val = { 1,2,0,3,4,0,5,6,0 };

	Matrix_Wrapper mw = { num_rows, num_columns, val.data(), leading_diemnsion };
	mw.transpose();

	constexpr auto ref_num_rows = 2;
	constexpr auto ref_num_columns = 3;
	std::vector<double> ref_val = { 1,3,5,2,4,6 };
	Matrix_Wrapper ref = { ref_num_rows, ref_num_columns, ref_val.data() };

	EXPECT_EQ(mw, ref);
}
