#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

using stan::math::check_finite;

TEST(ErrorHandlingScalar,CheckFinite) {
  const char* function = "check_finite";
  double x = 0;
 
  EXPECT_NO_THROW(check_finite(function, "x", x))
    << "check_finite should be true with finite x: " << x;
  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
    << "check_finite should throw exception on Inf: " << x;
  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error) 
    << "check_finite should throw exception on -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
    << "check_finite should throw exception on NaN: " << x;
}

TEST(ErrorHandlingScalar,CheckFinite_nan) {
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_finite(function, "x", nan), std::domain_error);
}
