#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  std::vector<var> x;
  x.push_back(var(1.0));
  x.push_back(var(2.0));
  x.push_back(var(3.0));

  EXPECT_NO_THROW(check_positive(function, "x", x));

  stan::math::recover_memory();
}
