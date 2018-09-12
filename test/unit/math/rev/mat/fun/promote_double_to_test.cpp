#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(AgradRevMatrix, promote_double_to_int) {
  int x = 3;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(AgradRevMatrix, promote_double_to_double) {
  double x = 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<0>(y))>::value));
  EXPECT_FLOAT_EQ(x, std::get<0>(y).val());
}

TEST(AgradRevMatrix, promote_double_to_mix) {
  auto y1 = stan::math::promote_double_to<stan::math::var>(
      std::make_tuple(1.0, 2, 3.0));
  auto y2 = stan::math::promote_double_to<stan::math::var>(
      std::make_tuple(1.0, 2));
  EXPECT_TRUE((std::is_same<std::tuple<stan::math::var, stan::math::var>,
                            decltype(y1)>::value));
  EXPECT_TRUE((std::is_same<std::tuple<stan::math::var>, decltype(y2)>::value));
  EXPECT_EQ(2, std::tuple_size<decltype(y1)>::value);
  EXPECT_EQ(1, std::tuple_size<decltype(y2)>::value);
  EXPECT_FLOAT_EQ(1.0, std::get<0>(y1).val());
  EXPECT_FLOAT_EQ(3.0, std::get<1>(y1).val());
  EXPECT_FLOAT_EQ(1.0, std::get<0>(y2).val());
}

TEST(AgradRevMatrix, promote_double_to_std_vector_int) {
  std::vector<int> x(3);
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));

  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(AgradRevMatrix, promote_double_to_std_vector_double) {
  std::vector<double> x = {{1.0, 2.0}};
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));

  EXPECT_TRUE((std::is_same<std::tuple<std::vector<stan::math::var> >,
                            decltype(y)>::value));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_FLOAT_EQ(1.0, std::get<0>(y)[0].val());
  EXPECT_FLOAT_EQ(2.0, std::get<0>(y)[1].val());
}

TEST(AgradRevMatrix, promote_double_to_std_vector_mix) {
  std::vector<int> xi = {{1, 2}};
  std::vector<double> xd1 = {{3.0, 4.0}};
  std::vector<double> xd2 = {{5.0, 6.0}};
  auto y1 = stan::math::promote_double_to<stan::math::var>(
      std::make_tuple(xd1, xi, xd2));
  auto y2 = stan::math::promote_double_to<stan::math::var>(
      std::make_tuple(xd1, xi));
  EXPECT_TRUE((std::is_same<std::tuple<std::vector<stan::math::var>,
                                       std::vector<stan::math::var> >,
                            decltype(y1)>::value));
  EXPECT_TRUE((std::is_same<std::tuple<std::vector<stan::math::var> >,
                            decltype(y2)>::value));
  EXPECT_EQ(2, std::tuple_size<decltype(y1)>::value);
  EXPECT_EQ(1, std::tuple_size<decltype(y2)>::value);
  EXPECT_FLOAT_EQ(3.0, std::get<0>(y1)[0].val());
  EXPECT_FLOAT_EQ(4.0, std::get<0>(y1)[1].val());
  EXPECT_FLOAT_EQ(5.0, std::get<1>(y1)[0].val());
  EXPECT_FLOAT_EQ(6.0, std::get<1>(y1)[1].val());
  EXPECT_FLOAT_EQ(3.0, std::get<0>(y2)[0].val());
  EXPECT_FLOAT_EQ(4.0, std::get<0>(y2)[1].val());
}

TEST(AgradRevMatrix, promote_double_to_VectorXf) {
  Eigen::VectorXf x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(AgradRevMatrix, promote_double_to_RowVectorXf) {
  Eigen::RowVectorXf x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}
TEST(AgradRevMatrix, promote_double_to_MatrixXf) {
  Eigen::MatrixXf x(1, 1);
  x << 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(0, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE((std::is_same<std::tuple<>, decltype(y)>::value));
}

TEST(AgradRevMatrix, promote_double_to_VectorXd) {
  Eigen::VectorXd x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_FLOAT_EQ(x(0), std::get<0>(y)(0).val());
}

TEST(AgradRevMatrix, promote_double_to_RowVectorXd) {
  Eigen::RowVectorXd x(1);
  x << 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_FLOAT_EQ(x(0), std::get<0>(y)(0).val());
}

TEST(AgradRevMatrix, promote_double_to_MatrixXd) {
  Eigen::MatrixXd x(1, 1);
  x << 3.0;
  auto y = stan::math::promote_double_to<stan::math::var>(std::make_tuple(x));
  EXPECT_EQ(1, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_FLOAT_EQ(x(0), std::get<0>(y)(0).val());
}

TEST(AgradRevMatrix, promote_double_to_Eigen_Mix) {
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

  auto y = stan::math::promote_double_to<stan::math::var>(
      std::make_tuple(x0, x1, x2, x3, x4, x5));
  EXPECT_EQ(3, std::tuple_size<decltype(y)>::value);
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<0>(y)(0))>::value));
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<1>(y)(0))>::value));
  EXPECT_TRUE(
      (std::is_same<stan::math::var&, decltype(std::get<2>(y)(0))>::value));
  EXPECT_FLOAT_EQ(x1(0), std::get<0>(y)(0).val());
  EXPECT_FLOAT_EQ(x3(0), std::get<1>(y)(0).val());
  EXPECT_FLOAT_EQ(x5(0), std::get<2>(y)(0).val());
}
