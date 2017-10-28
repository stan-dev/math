#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingScalar, CheckPositive) {
  using stan::math::check_positive;
  const std::string function = "check_positive";

  EXPECT_NO_THROW(check_positive(function, "x", 3.0));
}

TEST(ErrorHandlingScalar, CheckPositive_nan) {
  using stan::math::check_positive;
  const std::string function = "check_positive";

  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_positive(function, "x", nan),
               std::domain_error);
}
