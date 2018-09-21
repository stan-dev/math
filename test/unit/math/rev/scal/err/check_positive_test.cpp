#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <exception>
#include <limits>

using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  var x = std::numeric_limits<var>::quiet_NaN();
  EXPECT_THROW(check_positive(function, "x", x), std::domain_error);

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckPositiveVarCheckUnivariate) {
  using stan::math::check_positive;
  using stan::math::var;

  const char* function = "check_positive";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::instance().var_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_NO_THROW(check_positive(function, "a", a));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  stan::math::recover_memory();
}
