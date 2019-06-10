#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

using stan::math::check_positive_finite;

TEST(ErrorHandlingScalar, CheckPositiveFinite_Vector) {
  const char* function = "check_positive_finite";
  std::vector<double> x = {1.5, 0.1, 1};
  ASSERT_NO_THROW(check_positive_finite(function, "x", x))
      << "check_positive_finite should be true with finite x";

  x = {1, 2, std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on Inf";

  x = {-1, 2, std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on negative x";

  x = {0, 2, std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on x=0";

  x = {1, 2, -std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on -Inf";

  x = {1, 2, std::numeric_limits<double>::quiet_NaN()};
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on NaN";
}

TEST(ErrorHandlingScalar, CheckPositiveFinite_nan) {
  const char* function = "check_positive_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = nan;
    EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error);
    x[i] = i;
  }
}
