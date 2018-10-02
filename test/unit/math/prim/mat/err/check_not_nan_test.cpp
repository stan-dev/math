#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <limits>

TEST(ErrorHandlingMatrix, checkNotNanEigenRow) {
  stan::math::vector_d y;
  y.resize(3);

  EXPECT_NO_THROW(stan::math::check_not_nan("checkNotNanEigenRow(%1)", "y", y));
  EXPECT_NO_THROW(stan::math::check_not_nan("checkNotNanEigenRow(%1)", "y", y));

  y(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::check_not_nan("checkNotNanEigenRow", "y", y),
               std::domain_error);
  EXPECT_THROW(stan::math::check_not_nan("checkNotNanEigenRow", "y", y),
               std::domain_error);
}
