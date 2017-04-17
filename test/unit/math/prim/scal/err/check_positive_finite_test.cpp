#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

using stan::math::check_positive_finite;

TEST(ErrorHandlingScalar,CheckPositiveFinite) {
  const char* function = "check_positive_finite";
  double x = 1;
 
  EXPECT_NO_THROW(check_positive_finite(function, "x", x))
    << "check_positive_finite should be true with finite x: " << x;
  x = -1;
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on x= " << x;
  x = 0;
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on x= " << x;
  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on Inf: " << x;
  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error) 
    << "check_positive_finite should throw exception on -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
    << "check_positive_finite should throw exception on NaN: " << x;
}

TEST(ErrorHandlingScalar,CheckPositiveFinite_nan) {
  const char* function = "check_positive_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_positive_finite(function, "x", nan),
               std::domain_error);
}
