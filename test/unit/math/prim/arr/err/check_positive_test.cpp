#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandlingScalar, CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    EXPECT_NO_THROW(check_positive(function, "x", x));
  }
}

TEST(ErrorHandlingScalar, CheckPositive_nan) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = nan;
    EXPECT_THROW(check_positive(function, "x", x), std::domain_error);
    x[i] = i;
  }
}
