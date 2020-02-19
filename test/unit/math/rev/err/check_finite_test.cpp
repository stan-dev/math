#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(AgradRevErrorHandlingScalar, CheckFiniteVarCheckVectorized) {
  using stan::math::check_finite;
  using stan::math::var;
  using std::vector;

  int N = 5;
  const char* function = "check_finite";
  vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();

  EXPECT_EQ(5U, stack_size);
  EXPECT_NO_THROW(check_finite(function, "a", a));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(5U, stack_size_after_call);

  a[1] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "a", a), std::domain_error);
  stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(6U, stack_size_after_call);

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckFinite) {
  using stan::math::check_finite;
  using stan::math::var;

  const char* function = "check_bounded";
  const char* name = "x";
  var x = 0;

  EXPECT_NO_THROW(check_finite(function, name, x))
      << "check_finite should be TRUE with x: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_finite(function, name, x), std::domain_error)
      << "check_finite should throw with x: " << x;

  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, name, x), std::domain_error)
      << "check_finite should throw with x: " << x;

  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, name, x), std::domain_error)
      << "check_finite should throw with x: " << x;
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckFiniteVarCheckUnivariate) {
  using stan::math::check_finite;
  using stan::math::var;

  const char* function = "check_finite";
  var a(5.0);

  size_t stack_size
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_NO_THROW(check_finite(function, "a", a));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  a = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "a", a), std::domain_error);
  stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(2U, stack_size_after_call);

  stan::math::recover_memory();
}
