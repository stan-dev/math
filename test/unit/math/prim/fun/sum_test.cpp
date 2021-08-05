#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

TEST(MathFunctions, sumZeroSize) {
  std::vector<double> x;
  EXPECT_FLOAT_EQ(0.0, stan::math::sum(x));
}

TEST(MathFunctions, sum) {
  std::vector<double> x(3);

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;

  EXPECT_FLOAT_EQ(6.0, stan::math::sum(x));
}

TEST(MathFunctions, sum_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(3);

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = nan;

  EXPECT_TRUE(std::isnan(stan::math::sum(x)));
}
TEST(MathMatrixPrimArr, sum_vector_int) {
  std::vector<int> x(3);
  EXPECT_EQ(0, stan::math::sum(x));
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  EXPECT_EQ(6, stan::math::sum(x));
}

TEST(MathMatrixPrimMat, sumVector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::sum;

  stan::math::vector_d v;
  EXPECT_FLOAT_EQ(0.0, sum(v));

  v = stan::math::vector_d(1);
  v[0] = 5.0;
  EXPECT_FLOAT_EQ(5.0, sum(v));

  v = stan::math::vector_d(3);
  v[0] = 5.0;
  v[1] = 10.0;
  v[2] = 100.0;
  EXPECT_FLOAT_EQ(115.0, sum(v));
}

TEST(MathMatrixPrimMat, sumRowVector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::sum;

  stan::math::row_vector_d rv;
  EXPECT_FLOAT_EQ(0.0, sum(rv));

  rv = stan::math::row_vector_d(1);
  rv[0] = 5.0;
  EXPECT_FLOAT_EQ(5.0, sum(rv));

  rv = stan::math::row_vector_d(3);
  rv[0] = 5.0;
  rv[1] = 10.0;
  rv[2] = 100.0;
  EXPECT_FLOAT_EQ(115.0, sum(rv));
}

TEST(MathMatrixPrimMat, sumMatrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::sum;
  stan::math::matrix_d m;
  EXPECT_FLOAT_EQ(0.0, sum(m));

  m = stan::math::matrix_d(1, 1);
  m << 5.0;
  EXPECT_FLOAT_EQ(5.0, sum(m));

  m = stan::math::matrix_d(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_FLOAT_EQ(21.0, sum(m));
}

template <typename T>
using sum_return_t = decltype(stan::math::sum(std::declval<T>()));

TEST(MathMatrixPrimMat, sumIsTemplated) {
  using Eigen::Matrix;
  EXPECT_SAME_TYPE(int, sum_return_t<Matrix<int, 2, 3>>);
  EXPECT_SAME_TYPE(double, sum_return_t<Matrix<double, 4, 2>>);
  EXPECT_SAME_TYPE(float, sum_return_t<std::vector<float>>);
}
