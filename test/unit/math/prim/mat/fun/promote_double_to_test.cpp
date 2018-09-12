#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(MathFunctions, promote_double_to_int) {
  int x = 3;
  auto y = stan::math::promote_double_to<double>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(MathFunctions, promote_double_to_double) {
  double x = 3.0;
  auto y = stan::math::promote_double_to<float>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<0>(y))>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x), std::get<0>(y));
}

TEST(MathFunctions, promote_double_to_mix) {
  auto y1 = stan::math::promote_double_to<float>(std::make_tuple(1.0, 2, 3.0));
  auto y2 = stan::math::promote_double_to<float>(std::make_tuple(1.0, 2));
  EXPECT_TRUE((std::is_same<std::tuple<float, float>, decltype(y1)>::value));
  EXPECT_TRUE((std::is_same<std::tuple<float>, decltype(y2)>::value));
  EXPECT_EQ(2, std::tuple_size<decltype(y1)>::value);
  EXPECT_EQ(1, std::tuple_size<decltype(y2)>::value);
  EXPECT_FLOAT_EQ(static_cast<float>(1.0), std::get<0>(y1));
  EXPECT_FLOAT_EQ(static_cast<float>(3.0), std::get<1>(y1));
  EXPECT_FLOAT_EQ(static_cast<float>(1.0), std::get<0>(y2));
}

TEST(MathFunctions, promote_double_to_std_vector_int) {
  std::vector<int> x(3);
  auto y = stan::math::promote_double_to<double>(std::make_tuple(x));

  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(MathFunctions, promote_double_to_std_vector_double) {
  std::vector<double> x = {{1.0, 2.0}};
  auto y = stan::math::promote_double_to<float>(std::make_tuple(x));

  EXPECT_TRUE(
      (std::is_same<std::tuple<std::vector<float> >, decltype(y)>::value));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_FLOAT_EQ(static_cast<float>(1.0), std::get<0>(y)[0]);
  EXPECT_FLOAT_EQ(static_cast<float>(2.0), std::get<0>(y)[1]);
}

TEST(MathFunctions, promote_double_to_std_vector_mix) {
  std::vector<int> xi = {{1, 2}};
  std::vector<double> xd1 = {{3.0, 4.0}};
  std::vector<double> xd2 = {{5.0, 6.0}};
  auto y1 = stan::math::promote_double_to<float>(std::make_tuple(xd1, xi, xd2));
  auto y2 = stan::math::promote_double_to<float>(std::make_tuple(xd1, xi));
  EXPECT_TRUE((std::is_same<std::tuple<std::vector<float>, std::vector<float> >,
                            decltype(y1)>::value));
  EXPECT_TRUE(
      (std::is_same<std::tuple<std::vector<float> >, decltype(y2)>::value));
  EXPECT_EQ(2, std::tuple_size<decltype(y1)>::value);
  EXPECT_EQ(1, std::tuple_size<decltype(y2)>::value);
  EXPECT_FLOAT_EQ(static_cast<float>(3.0), std::get<0>(y1)[0]);
  EXPECT_FLOAT_EQ(static_cast<float>(4.0), std::get<0>(y1)[1]);
  EXPECT_FLOAT_EQ(static_cast<float>(5.0), std::get<1>(y1)[0]);
  EXPECT_FLOAT_EQ(static_cast<float>(6.0), std::get<1>(y1)[1]);
  EXPECT_FLOAT_EQ(static_cast<float>(3.0), std::get<0>(y2)[0]);
  EXPECT_FLOAT_EQ(static_cast<float>(4.0), std::get<0>(y2)[1]);
}

TEST(MathFunctions, promote_double_to_VectorXf) {
  Eigen::VectorXf x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<double>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(MathFunctions, promote_double_to_RowVectorXf) {
  Eigen::RowVectorXf x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<double>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}
TEST(MathFunctions, promote_double_to_MatrixXf) {
  Eigen::MatrixXf x(1, 1);
  x << 3.0;
  auto y = stan::math::promote_double_to<double>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(MathFunctions, promote_double_to_VectorXd) {
  Eigen::VectorXd x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<float>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x(0)), std::get<0>(y)(0));
}

TEST(MathFunctions, promote_double_to_RowVectorXd) {
  Eigen::RowVectorXd x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<float>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x(0)), std::get<0>(y)(0));
}

TEST(MathFunctions, promote_double_to_MatrixXd) {
  Eigen::MatrixXd x(1, 1);
  x << 3.0;
  auto y = stan::math::promote_double_to<float>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x(0)), std::get<0>(y)(0));
}

TEST(MathFunctions, promote_double_to_Eigen_Mix) {
  Eigen::VectorXf x0(1);
  Eigen::RowVectorXf x2(1);
  Eigen::MatrixXf x4(1, 1);
  Eigen::VectorXd x1(1);
  Eigen::RowVectorXd x3(1);
  Eigen::MatrixXd x5(1, 1);
  x0 << 0.0;
  x1 << 1.0;
  x2 << 2.0;
  x3 << 3.0;
  x4 << 4.0;
  x5 << 5.0;

  auto y = stan::math::promote_double_to<float>(
      std::make_tuple(x0, x1, x2, x3, x4, x5));
  EXPECT_EQ(3, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<1>(y)(0))>::value));
  EXPECT_TRUE((std::is_same<float&, decltype(std::get<2>(y)(0))>::value));
  EXPECT_FLOAT_EQ(static_cast<float>(x1(0)), std::get<0>(y)(0));
  EXPECT_FLOAT_EQ(static_cast<float>(x3(0)), std::get<1>(y)(0));
  EXPECT_FLOAT_EQ(static_cast<float>(x5(0)), std::get<2>(y)(0));
}
