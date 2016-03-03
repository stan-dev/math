#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

using stan::math::check_equal;

TEST(ErrorHandlingScalar,CheckEqual) {
  const char* function = "check_equal";
  double x = 0.0;
  double eq = 0.0;
 
  EXPECT_TRUE(check_equal(function, "x", x, eq))
    << "check_equal should be true with x = eq";
  
  x = -1.0;
  EXPECT_THROW(check_equal(function, "x", x, eq),
               std::domain_error)
    << "check_equal should throw an exception with x < eq";

  x = eq;
  EXPECT_NO_THROW(check_equal(function, "x", x, eq))
    << "check_equal should not throw an exception with x == eq";

  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function, "x", x, eq), 
               std::domain_error)
    << "check_equal should be false with x == Inf and eq = 0.0";

  x = 10.0;
  eq = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function, "x", x, eq),
               std::domain_error)
    << "check_equal should throw an exception with x == 10.0 and eq == Inf";

  x = std::numeric_limits<double>::infinity();
  eq = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_equal(function, "x", x, eq))
    << "check_equal should not throw an exception with x == Inf and eq == Inf";
}


TEST(ErrorHandlingScalar,CheckEqual_nan) {
  const char* function = "check_equal";
  double x = 0.0;
  double eq = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_equal(function, "x", x, nan),
               std::domain_error);
  EXPECT_THROW(check_equal(function, "x", nan, eq),
               std::domain_error);
  EXPECT_THROW(check_equal(function, "x", nan, nan),
               std::domain_error);
}
