#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

using stan::math::var;
using stan::math::check_greater_or_equal;

TEST(AgradRevErrorHandlingScalar,CheckGreaterOrEqual) {
  const char* function = "check_greater_or_equal";
  var x = 10.0;
  var lb = 0.0;
 
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb)) 
    << "check_greater_or_equal should be true with x > lb";
  
  x = -1.0;
  EXPECT_THROW(check_greater_or_equal(function, "x", x, lb),
               std::domain_error)
    << "check_greater_or_equal should throw an exception with x < lb";

  x = lb;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
    << "check_greater_or_equal should not throw an exception with x == lb";

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
    << "check_greater should be true with x == Inf and lb = 0.0";

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater_or_equal(function, "x", x, lb),
               std::domain_error)
    << "check_greater should throw an exception with x == 10.0 and lb == Inf";

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
    << "check_greater should not throw an exception with x == Inf and lb == Inf";
  stan::math::recover_memory();
}


TEST(AgradRevErrorHandlingScalar, CheckGreaterOrEqualVarCheckUnivariate) {
  using stan::math::var;
  using stan::math::check_greater_or_equal;

  const char* function = "check_greater_or_equal";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(1U,stack_size);
  EXPECT_NO_THROW(check_greater_or_equal(function,"a",a,2.0));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(1U,stack_size_after_call);

  EXPECT_THROW(check_greater_or_equal(function,"a",a,10.0),std::domain_error);
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(1U,stack_size_after_call);

  stan::math::recover_memory();
}
