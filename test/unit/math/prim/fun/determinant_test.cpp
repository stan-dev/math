#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, dimensionValidation) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::determinant;
  Matrix<double, Dynamic, Dynamic> x(3, 3);
  x << 1, 2, 3, 1, 4, 9, 1, 8, 27;

  ASSERT_FALSE(stan::math::is_nan(determinant(x)));

  Matrix<double, Dynamic, Dynamic> xx(3, 2);
  xx << 1, 2, 3, 1, 4, 9;
  EXPECT_THROW(stan::math::determinant(xx), std::invalid_argument);
}
