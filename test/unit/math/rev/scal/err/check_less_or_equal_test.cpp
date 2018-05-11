#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_less_or_equal;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckLessOrEqual) {
  const char* function = "check_less_or_equal";
  var x = -10.0;
  var lb = 0.0;

  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less_or_equal should be true with x < lb";

  x = 1.0;
  EXPECT_THROW(check_less_or_equal(function, "x", x, lb), std::domain_error)
      << "check_less_or_equal should throw an exception with x > lb";

  x = lb;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less_or_equal should not throw an exception with x == lb";

  x = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less should be true with x == -Inf and lb = 0.0";

  x = -10.0;
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less_or_equal(function, "x", x, lb), std::domain_error)
      << "check_less should throw an exception with x == -10.0 and lb == -Inf";

  x = -std::numeric_limits<double>::infinity();
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less should not throw an exception with x == -Inf and lb == "
         "-Inf";
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckLessOrEqualVarCheckUnivariate) {
  using stan::math::check_less_or_equal;
  using stan::math::var;

  const char* function = "check_less_or_equal";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::instance().var_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_THROW(check_less_or_equal(function, "a", a, 2.0), std::domain_error);

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  EXPECT_NO_THROW(check_less_or_equal(function, "a", a, 5.0));

  stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  EXPECT_NO_THROW(check_less_or_equal(function, "a", a, 10.0));
  stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  stan::math::recover_memory();
}
