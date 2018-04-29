#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRevErrorHandlingScalar, CheckConsistentSizeVarCheckVectorized) {
  using stan::math::check_consistent_sizes;
  using stan::math::var;
  using std::vector;

  int N = 5;
  const char* function = "check_consistent_size";
  vector<var> a;
  vector<var> b;

  for (int i = 0; i < N; ++i) {
    b.push_back(var(i + 1));
    a.push_back(var(i));
  }

  size_t stack_size = stan::math::ChainableStack::instance().var_stack_.size();

  EXPECT_EQ(10U, stack_size);
  EXPECT_NO_THROW(check_consistent_sizes(function, "a", a, "b", b));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance().var_stack_.size();
  EXPECT_EQ(10U, stack_size_after_call);
  stan::math::recover_memory();
}
