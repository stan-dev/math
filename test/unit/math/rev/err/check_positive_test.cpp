#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <exception>
#include <limits>
#include <vector>

TEST(AgradRevErrorHandlingArray, CheckPositive) {
  using stan::math::check_positive;
  using stan::math::var;

  const char* function = "check_positive";

  std::vector<var> x;
  x.push_back(var(1.0));
  x.push_back(var(2.0));
  x.push_back(var(3.0));

  EXPECT_NO_THROW(check_positive(function, "x", x));

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingMatrix, CheckPositive) {
  using stan::math::check_positive;
  using stan::math::var;

  const char* function = "check_positive";

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_mat(3);
  x_mat << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    EXPECT_NO_THROW(check_positive(function, "x", x_mat));
  }

  x_mat(0) = 0;

  EXPECT_THROW(check_positive(function, "x", x_mat), std::domain_error);

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckPositive) {
  using stan::math::check_positive;
  using stan::math::var;

  const char* function = "check_positive";

  var x = std::numeric_limits<var>::quiet_NaN();
  EXPECT_THROW(check_positive(function, "x", x), std::domain_error);

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckPositiveVarCheckUnivariate) {
  using stan::math::check_positive;
  using stan::math::var;

  const char* function = "check_positive";
  var a(5.0);

  size_t stack_size
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_NO_THROW(check_positive(function, "a", a));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  stan::math::recover_memory();
}
