#include <gtest/gtest.h>
#include <stan/math/rev/arr.hpp>
#include <string>
#include <vector>

TEST(AgradRevErrorHandlingScalar, CheckNotNanVarCheckVectorized) {
  using stan::math::var;
  using std::vector;
  using stan::math::check_not_nan;

  int N = 5;
  const std::string function = "check_not_nan";
  vector<var> a;

  for (int i = 0; i < N; ++i) a.push_back(var(i));

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(5U, stack_size);
  EXPECT_NO_THROW(check_not_nan(function, "a", a));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U, stack_size_after_call);
  stan::math::recover_memory();
}

TEST(ErrorHandlingScalar, CheckNotNanVarCheckVectorized) {
  using stan::math::var;
  using std::vector;
  using stan::math::check_not_nan;

  int N = 5;
  const std::string function = "check_not_nan";
  vector<var> a;

  for (int i = 0; i < N; ++i) a.push_back(var(i));

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(5U, stack_size);
  EXPECT_NO_THROW(check_not_nan(function, "a", a));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U, stack_size_after_call);
  stan::math::recover_memory();
}
