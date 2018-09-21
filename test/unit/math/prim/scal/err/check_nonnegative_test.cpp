#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_nonnegative;

TEST(ErrorHandlingScalar, CheckNonnegative) {
  const char* function = "check_nonnegative";
  double x = 0;

  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative should be true with finite x: " << x;

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative should be true with x = Inf: " << x;

  x = -0.01;
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = " << x;

  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception on NaN: " << x;
}

TEST(ErrorHandlingScalar, CheckNonnegative_nan) {
  const char* function = "check_nonnegative";
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_nonnegative(function, "x", nan), std::domain_error);
}
