#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(AgradRevErrorHandlingScalar, CheckLessOrEqualVarCheckVectorized) {
  using stan::math::check_less_or_equal;
  using stan::math::var;

  int N = 5;
  const char* function = "check_less_or_equal";
  std::vector<var> a;

  for (int i = 0; i < N; ++i)
    a.push_back(var(i));

  size_t stack_size
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();

  EXPECT_EQ(5U, stack_size);
  EXPECT_NO_THROW(check_less_or_equal(function, "a", a, 10.0));

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(5U, stack_size_after_call);

  EXPECT_THROW(check_less_or_equal(function, "a", a, 2.0), std::domain_error);
  stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(5U, stack_size_after_call);

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingMatrix, CheckLessOrEqual_Matrix) {
  using stan::math::check_less_or_equal;
  using stan::math::var;
  const char* function = "check_less_or_equal";
  var x;
  var high;
  Eigen::Matrix<var, Eigen::Dynamic, 1> x_vec;
  Eigen::Matrix<var, Eigen::Dynamic, 1> high_vec;
  x_vec.resize(3);
  high_vec.resize(3);

  // x_vec, high
  x_vec << -5, 0, 5;
  high = 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  x_vec << -5, 0, 5;
  high = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  x_vec << -5, 0, 5;
  high = 5;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high = 5;
  EXPECT_THROW(check_less_or_equal(function, "x", x_vec, high),
               std::domain_error);

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  // x_vec, high_vec
  x_vec << -5, 0, 5;
  high_vec << 0, 5, 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  x_vec << -5, 0, 5;
  high_vec << std::numeric_limits<double>::infinity(), 10, 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  x_vec << -5, 0, 5;
  high_vec << 10, 10, 5;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high_vec << 10, 10, 10;
  EXPECT_THROW(check_less_or_equal(function, "x", x_vec, high_vec),
               std::domain_error);

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high_vec << 10, 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  // x, high_vec
  x = -100;
  high_vec << 0, 5, 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));

  x = 10;
  high_vec << 100, 200, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));

  x = 5;
  high_vec << 100, 200, 5;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));

  x = std::numeric_limits<double>::infinity();
  high_vec << 10, 20, 30;
  EXPECT_THROW(check_less_or_equal(function, "x", x, high_vec),
               std::domain_error);

  x = std::numeric_limits<double>::infinity();
  high_vec << std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckLessOrEqual) {
  using stan::math::check_less_or_equal;
  using stan::math::var;

  const char* function = "check_less_or_equal";
  var x = -10.0;
  var lb = 0.0;

  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less_or_equal should be true with x < lb";

  x = 1.0;
  EXPECT_THROW(check_less_or_equal(function, "x", x, lb), std::domain_error)
      << "check_less_or_equal should throw an exception with x > lb";

  x = lb;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less_or_equal should not throw an exception with x == lb";

  x = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less should be true with x == -Inf and lb = 0.0";

  x = -10.0;
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less_or_equal(function, "x", x, lb), std::domain_error)
      << "check_less should throw an exception with x == -10.0 and lb == -Inf";

  x = -std::numeric_limits<double>::infinity();
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, lb))
      << "check_less should not throw an exception with x == -Inf and lb == "
         "-Inf";
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckLessOrEqualVarCheckUnivariate) {
  using stan::math::check_less_or_equal;
  using stan::math::var;

  const char* function = "check_less_or_equal";
  var a(5.0);

  size_t stack_size
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();

  EXPECT_EQ(1U, stack_size);
  EXPECT_THROW(check_less_or_equal(function, "a", a, 2.0), std::domain_error);

  size_t stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  EXPECT_NO_THROW(check_less_or_equal(function, "a", a, 5.0));

  stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  EXPECT_NO_THROW(check_less_or_equal(function, "a", a, 10.0));
  stack_size_after_call
      = stan::math::ChainableStack::instance_->var_stack_.size()
        + stan::math::ChainableStack::instance_->var_nochain_stack_.size();
  EXPECT_EQ(1U, stack_size_after_call);

  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingMatrix, CheckLessOrEqualStdVecVarMatrix) {
  using stan::math::var;
  std::vector<std::vector<
      stan::conditional_var_value_t<var, Eigen::Matrix<var, -1, -1>>>>
      ar_mat = std::vector<std::vector<
          stan::conditional_var_value_t<var, Eigen::Matrix<var, -1, -1>>>>(
          4,
          std::vector<
              stan::conditional_var_value_t<var, Eigen::Matrix<var, -1, -1>>>(
              5, stan::conditional_var_value_t<var, Eigen::Matrix<var, -1, -1>>(
                     Eigen::Matrix<double, -1, -1>::Constant(
                         2, 3, std::numeric_limits<double>::quiet_NaN()))));
  EXPECT_THROW(stan::math::check_less_or_equal("test_g", "ar_mat", ar_mat, 0),
               std::domain_error);
}
