#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_nonnegative;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckNonnegative) {
  const char* function = "check_nonnegative";
  var x = 0;

  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative should be true with finite x: " << x;

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative should be true with x = Inf: " << x;

  x = -0.01;
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = " << x;

  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception on NaN: " << x;
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckNonnegativeVarCheckUnivariate) {
  using stan::math::check_nonnegative;
  using stan::math::var;

  const char* function = "check_nonnegative";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::instance().var_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_NO_THROW(check_nonnegative(function, "a", a));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  a = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_nonnegative(function, "a", a));
  stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(2U, stack_size_after_call);

  a = 0.0;
  EXPECT_NO_THROW(check_nonnegative(function, "a", a));
  stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(3U, stack_size_after_call);

  a = -1.1;
  EXPECT_THROW(check_nonnegative(function, "a", a), std::domain_error);
  stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(4U, stack_size_after_call);
  stan::math::recover_memory();
}
