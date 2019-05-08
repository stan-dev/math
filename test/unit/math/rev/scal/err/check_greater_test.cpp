#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_greater;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckGreater) {
  const char* function = "check_greater";
  var x = 10.0;
  var lb = 0.0;

  EXPECT_NO_THROW(check_greater(function, "x", x, lb))
      << "check_greater should be true with x > lb";

  x = -1.0;
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x < lb";

  x = lb;
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == lb";

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater(function, "x", x, lb))
      << "check_greater should be true with x == Inf and lb = 0.0";

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == 10.0 and lb == Inf";

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == Inf and lb == Inf";
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckGreaterVarCheckUnivariate) {
  using stan::math::check_greater;
  using stan::math::var;

  const char* function = "check_greater";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::instance().var_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_NO_THROW(check_greater(function, "a", a, 2.0));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  EXPECT_THROW(check_greater(function, "a", a, 10.0), std::domain_error);
  stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  stan::math::recover_memory();
}
