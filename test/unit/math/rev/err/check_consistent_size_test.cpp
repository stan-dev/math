#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(AgradRevErrorHandlingScalar, CheckConsistentSizeVarCheckVectorized) {
  using stan::math::check_consistent_size;
  using stan::math::var;
  using std::vector;

  int N = 5;
  const char* function = "check_consistent_size";
  vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();

  EXPECT_EQ(5U, stack_size);
  EXPECT_NO_THROW(check_consistent_size(function, "a", a, 5U));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(5U, stack_size_after_call);
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingMatrix, checkConsistentSize) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::check_consistent_size;
  using stan::math::size;
  using stan::math::var;

  const char* function = "check_consistent_size";
  const char* name1 = "name1";

  Matrix<var, Dynamic, 1> v1(4);
  v1 << 4.0, 5.0, 6.0, 7.0;
  EXPECT_EQ(4U, stan::math::size(v1));
  EXPECT_NO_THROW(check_consistent_size(function, name1, v1, 4U));
  EXPECT_THROW(check_consistent_size(function, name1, v1, 2U),
               std::invalid_argument);
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingMatrix, checkConsistentSize_nan) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::check_consistent_size;
  using stan::math::size;
  using stan::math::var;

  const char* function = "check_consistent_size";
  const char* name1 = "name1";

  double nan = std::numeric_limits<double>::quiet_NaN();

  Matrix<var, Dynamic, 1> v1(4);
  v1 << nan, nan, 4, nan;
  EXPECT_EQ(4U, stan::math::size(v1));
  EXPECT_NO_THROW(check_consistent_size(function, name1, v1, 4U));
  EXPECT_THROW(check_consistent_size(function, name1, v1, 2U),
               std::invalid_argument);
  stan::math::recover_memory();
}
