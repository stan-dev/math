#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradRevErrorHandlingScalar,CheckNotNan) {
  using stan::math::var;
  using stan::math::check_not_nan;
  const char* function = "check_not_nan";

  var x = 0;
  double x_d = 0;
 
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
    << "check_not_nan should be true with finite x: " << x;
  EXPECT_NO_THROW(check_not_nan(function, "x", x_d))
    << "check_not_nan should be true with finite x: " << x_d;
  
  x = std::numeric_limits<var>::infinity();
  x_d = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
    << "check_not_nan should be true with x = Inf: " << x;
  EXPECT_NO_THROW(check_not_nan(function, "x", x_d))
    << "check_not_nan should be true with x = Inf: " << x_d;

  x = -std::numeric_limits<var>::infinity();
  x_d = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
    << "check_not_nan should be true with x = -Inf: " << x;
  EXPECT_NO_THROW(check_not_nan(function, "x", x_d))
    << "check_not_nan should be true with x = -Inf: " << x_d;

  x = std::numeric_limits<var>::quiet_NaN();
  x_d = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_not_nan(function, "x", x), std::domain_error)
    << "check_not_nan should throw exception on NaN: " << x;
  EXPECT_THROW(check_not_nan(function, "x", x_d), std::domain_error)
    << "check_not_nan should throw exception on NaN: " << x_d;
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckNotNanVarCheckUnivariate) {
  using stan::math::var;
  using stan::math::check_not_nan;

  const char* function = "check_not_nan";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(1U,stack_size);
  EXPECT_NO_THROW(check_not_nan(function,"a",a));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(1U,stack_size_after_call);

  stan::math::recover_memory();
}

TEST(ErrorHandlingScalar, CheckNotNanVarCheckUnivariate) {
  using stan::math::var;
  using stan::math::check_not_nan;
  
  const char* function = "check_not_nan";
  var a(5.0);

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_TRUE(1U == stack_size);
  EXPECT_NO_THROW(check_not_nan(function,"a",a));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_TRUE(1U == stack_size_after_call);

  stan::math::recover_memory();
}

