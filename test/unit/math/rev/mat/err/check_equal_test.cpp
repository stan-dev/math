#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_equal;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar,CheckEqualMatrix) {
  const char* function = "check_equal";
  Eigen::Matrix<var,Eigen::Dynamic,1> x_vec;
  Eigen::Matrix<var,Eigen::Dynamic,1> eq_vec;
  x_vec.resize(3);
  eq_vec.resize(3);

  // x_vec, low_vec
  x_vec   << -1, 0, 1;
  eq_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_equal(function, "x", x_vec, eq_vec)) 
    << "check_equal: matrix<3,1>, matrix<3,1>";

  x_vec   <<   -1,    0,   1;
  eq_vec << -1.1, -0.1, 0.9;
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec),
               std::domain_error) 
    << "check_equal: matrix<3,1>, matrix<3,1>";
  
  x_vec   << -1, 0,  1;
  eq_vec << -2, -1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec), 
               std::domain_error) 
    << "check_equal: matrix<3,1>, matrix<3,1>, should fail with infinity";
  stan::math::recover_memory();
}

TEST(AgradRevErrorHandlingScalar, CheckEqualVarCheckMatrix) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  int N = 2;
  const char* function = "check_not_nan";
  Matrix<var,Dynamic,Dynamic> a(N,N);
  Matrix<var,Dynamic,Dynamic> b(N,N);

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j){
      a(i,j) = var(i+j+1);
      b(i,j) = var(i+j+1);
    }

  size_t stack_size = stan::math::ChainableStack::var_stack_.size();
  size_t stack_size_expected = 2 * N * N;

  EXPECT_EQ(stack_size_expected,stack_size);
  EXPECT_NO_THROW(check_equal(function,"a",a,b));

  size_t stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(stack_size_expected,stack_size_after_call);

  b(1,1) = 45;
  EXPECT_THROW(check_equal(function,"a",a,b),std::domain_error);
  stack_size_after_call = stan::math::ChainableStack::var_stack_.size();
  EXPECT_EQ(stack_size_expected + 1,stack_size_after_call);

  stan::math::recover_memory();
}
