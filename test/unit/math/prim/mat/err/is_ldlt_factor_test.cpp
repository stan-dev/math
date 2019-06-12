#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isLDLTFactor_nan) {
  using stan::math::is_ldlt_factor;

  double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x(2, 2);
  stan::math::LDLT_factor<double, -1, -1> ldlt_x;

  x << nan, 1, 1, 3;
  ldlt_x.compute(x);
  EXPECT_FALSE(ldlt_x.success());
  EXPECT_FALSE(is_ldlt_factor(ldlt_x));

  x << 3, nan, 1, 3;
  ldlt_x.compute(x);
  EXPECT_TRUE(ldlt_x.success());
  EXPECT_TRUE(is_ldlt_factor(ldlt_x));

  x << 3, 1, nan, 3;
  ldlt_x.compute(x);
  EXPECT_FALSE(ldlt_x.success());
  EXPECT_FALSE(is_ldlt_factor(ldlt_x));

  x << 3, 1, 1, nan;
  ldlt_x.compute(x);
  EXPECT_FALSE(ldlt_x.success());
  EXPECT_FALSE(is_ldlt_factor(ldlt_x));
}
