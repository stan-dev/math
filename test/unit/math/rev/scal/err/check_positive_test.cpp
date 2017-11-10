#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

using stan::math::var;

TEST(AgradRevErrorHandlingScalar,CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  EXPECT_NO_THROW(check_positive(function, "x", nan));

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckPositiveVarCheckUnivariate) {
  using stan::math::var;
  using stan::math::check_positive;

  const char* function = "check_positive";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(1U,stack_size);
  EXPECT_NO_THROW(check_positive(function,"a",a));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(1U,stack_size_after_call);

  stan::math::recover_memory();
}
