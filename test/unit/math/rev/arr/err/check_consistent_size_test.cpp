#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(AgradRevErrorHandlingScalar, CheckConsistentSizeVarCheckVectorized) {
  using stan::math::var;
  using std::vector;
  using stan::math::check_consistent_size;

  int N = 5;
  const char* function = "check_consistent_size";
  vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(5U,stack_size);
  EXPECT_NO_THROW(check_consistent_size(function,"a",a,5U));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(5U,stack_size_after_call);
  stan::math::recover_memory();
}
