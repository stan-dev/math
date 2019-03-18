#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

using stan::math::check_finite;

// ---------- check_finite: vector tests ----------
TEST(ErrorHandlingScalar, CheckFinite_Vector) {
  const char* function = "check_finite";
  std::vector<double> x = {-1, 0, 1};
  ASSERT_NO_THROW(check_finite(function, "x", x))
      << "check_finite should be true with finite x";

  x = {-1, 0, std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on Inf";

  x = {-1, 0, -std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on -Inf";

  x = {-1, 0, std::numeric_limits<double>::quiet_NaN()};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on NaN";
}

TEST(ErrorHandlingScalar, CheckFinite_nan) {
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {nan, 0, 1};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x = {1, nan, 1};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x = {1, 0, nan};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);
}
