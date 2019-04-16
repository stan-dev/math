#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, isNotNanEigenRow) {
  stan::math::vector_d y;
  y.resize(3);

  EXPECT_TRUE(stan::math::is_not_nan(y));
  EXPECT_TRUE(stan::math::is_not_nan(y));

  y(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(stan::math::is_not_nan(y));
  EXPECT_FALSE(stan::math::is_not_nan(y));
}
