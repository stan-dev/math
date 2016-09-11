#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(AgradRevErrorHandlingScalar, CheckBoundedVarCheckVectorized) {
  using stan::math::var;
  using std::vector;
  using stan::math::check_bounded;

  int N = 5;
  const char* function = "check_bounded";
  vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(5U,stack_size);
  EXPECT_NO_THROW(check_bounded(function,"a",a,-1.0,6.0));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U,stack_size_after_call);
  stan::math::recover_memory();
}
