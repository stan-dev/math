#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

using stan::math::check_greater;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckGreaterVarCheckVectorized) {
  using stan::math::var;
  using std::vector;
  using stan::math::check_greater;

  int N = 5;
  const char* function = "check_greater";
  vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(5U,stack_size);
  EXPECT_NO_THROW(check_greater(function,"a",a,-1.0));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U,stack_size_after_call);

  EXPECT_THROW(check_greater(function,"a",a,2.0),std::domain_error);
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U,stack_size_after_call);

  stan::math::recover_memory();
}
