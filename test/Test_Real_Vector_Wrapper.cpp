#include "gtest\gtest.h"
#include "../MS_Math_Library/src/Real_Vector.h"


//for google test cout message 
std::ostream& operator<<(std::ostream& os, const Real_Vector_Wrapper& vw)
{
	return os << vw.to_string();
}
//for google test cout message 

TEST(Real_Vector_Wrapper, construct_1)
{
	std::vector<double> vec1 = { 0,0,0 };
	
	Real_Vector_Wrapper vw1 = vec1;
	auto vw2 = vw1;

	EXPECT_EQ(vw1, vw2);
}
TEST(Real_Vector_Wrapper, operator_assign_1)
{
	std::vector<double> vec1 = { 0,0,0 };
	const std::vector<double> vec2 = { 1,2,3 };
	const std::vector<double> vec3 = { 6,5,4 };

	Real_Vector_Wrapper vw1 = vec1;
	Const_Real_Vector_Wrapper cvw2 = vec2;

	vw1.change_value(cvw2 + vec3);

	const std::vector<double> ref = { 7,7,7 };
	EXPECT_EQ(vec1, ref);
}
TEST(Real_Vector_Wrapper, operator_addition_assign_1)
{
	std::vector<double> vec1 = { 1,2,3 };
	const std::vector<double> vec2 = { 6,5,4 };

	Real_Vector_Wrapper vw1 = vec1;
	vw1 += vec2;

	const std::vector<double> ref = { 7,7,7 };
	EXPECT_EQ(vec1, ref);
}
TEST(Real_Vector_Wrapper, operator_addition_assign_2)
{
	std::vector<double> vec1 = { 1,2,3 };
	const std::vector<double> vec2 = { 6,5,4,3,2,1 };

	Real_Vector_Wrapper v1 = vec1;
	Const_Real_Vector_Wrapper cv2(vec2.data(), 3);
	v1 += cv2;

	const Real_Vector ref = { 7,7,7 };
	EXPECT_EQ(v1, ref);
}
TEST(Real_Vector_Wrapper, operator_scalar_multiplication_assign_1)
{
	std::vector<double> vec = { 1,2,3 };

	Real_Vector_Wrapper v = vec;
	v *= 2;

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(vec, ref);
}
TEST(Real_Vector_Wrapper, operator_scalar_multiplication_assign_2)
{
	std::vector<double> vec = { 1,2,3,4,5,6,7,8,9 };

	constexpr auto inc = 3;
	Real_Vector_Wrapper v = { vec, inc };
	v *= 2;

	const std::vector<double> ref = { 2,2,3,8,5,6,14,8,9 };
	EXPECT_EQ(vec, ref);
}
TEST(Real_Vector_Wrapper, inner_product_1)
{
	std::vector<double> vec1 = { 1,2,3 };
	const std::vector<double> vec2 = { 4,5,6 };

	const Real_Vector_Wrapper vw1 = vec1;
	const auto result = vw1.inner_product(vec2);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}

//Caution!
#ifdef _DEBUG
TEST(Real_Vector_Wrapper, different_size_1)
{
	std::vector<double> vec1 = { 1,2,3,4 };
	const std::vector<double> vec2 = { 1,2 };

	Real_Vector_Wrapper vw1 = vec1;
	EXPECT_ANY_THROW(vw1 += vec2);
}
#endif

TEST(Real_Vector_Wrapper, changed_vector_1)
{
	std::vector<double> vec1 = { 1,2,3 };

	Real_Vector_Wrapper cvw1 = vec1;
	vec1 = { 4,3,2 };

	EXPECT_EQ(cvw1, vec1);
}
TEST(Real_Vector_Wrapper, moved_vector_1)
{
	std::vector<double> vec = { 1,2,3 };
	Real_Vector_Wrapper v1 = vec;

	std::vector<double> new_vec = { 4,5,6,7,8,9 };
	vec = std::move(new_vec);//move may invoke changing data ptr

	const std::vector<double> ref = { 4,5,6,7,8,9 };
	EXPECT_NE(v1, ref);
}