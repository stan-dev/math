#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

using stan::math::check_equal;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckNotNanVarCheckVectorized) {
  using std::vector;

  int N = 5;
  const char* function = "check_not_nan";
  vector<var> a;
  vector<var> b;

  for (int i = 0; i < N; ++i){
   a.push_back(var(i));
   b.push_back(var(i));
  }

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();

  EXPECT_EQ(10U,stack_size);
  EXPECT_TRUE(check_equal(function,"a",a,b));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(10U,stack_size_after_call);

  b[1] = 4.0;
  EXPECT_THROW(check_equal(function,"a",a,b),std::domain_error);
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(11U,stack_size_after_call);

  stan::math::recover_memory();
}
