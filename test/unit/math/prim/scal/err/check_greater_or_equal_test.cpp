#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

using stan::math::check_greater_or_equal;

TEST(ErrorHandlingScalar,CheckGreaterOrEqual) {
  const char* function = "check_greater_or_equal";
  double x = 10.0;
  double lb = 0.0;
 
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb)) 
    << "check_greater_or_equal should be true with x > lb";
  
  x = -1.0;
  EXPECT_THROW(check_greater_or_equal(function, "x", x, lb),
               std::domain_error)
    << "check_greater_or_equal should throw an exception with x < lb";

  x = lb;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
    << "check_greater_or_equal should not throw an exception with x == lb";

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
    << "check_greater should be true with x == Inf and lb = 0.0";

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater_or_equal(function, "x", x, lb),
               std::domain_error)
    << "check_greater should throw an exception with x == 10.0 and lb == Inf";

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
    << "check_greater should not throw an exception with x == Inf and lb == Inf";
}

TEST(ErrorHandlingScalar,CheckGreaterOrEqual_nan) {
  const char* function = "check_greater_or_equal";
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_greater_or_equal(function, "x", nan, lb),
               std::domain_error);
  EXPECT_THROW(check_greater_or_equal(function, "x", x, nan),
               std::domain_error);
  EXPECT_THROW(check_greater_or_equal(function, "x", nan, nan),
               std::domain_error);
}
 
