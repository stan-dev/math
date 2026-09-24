#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <utility>

using stan::math::matrix_cl;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;

TEST(ScalarCl, kernel_expression_traits) {
  using stan::is_kernel_expression;
  using stan::is_kernel_expression_and_not_scalar;
  using stan::is_kernel_expression_lhs;
  using stan::is_matrix_cl;
  using stan::is_nonscalar_prim_or_rev_kernel_expression;
  using stan::is_prim_or_rev_kernel_expression;
  using stan::is_rev_kernel_expression;
  // A device scalar is a scalar operand of kernel generator expressions ...
  EXPECT_TRUE((is_kernel_expression<ScalarCl<double>>::value));
  EXPECT_TRUE((is_kernel_expression<const ScalarCl<double>&>::value));
  EXPECT_TRUE((is_prim_or_rev_kernel_expression<ScalarCl<double>>::value));
  // ... but never a matrix or an assignable matrix expression.
  EXPECT_FALSE((is_matrix_cl<ScalarCl<double>>::value));
  EXPECT_FALSE((is_kernel_expression_and_not_scalar<ScalarCl<double>>::value));
  EXPECT_FALSE((is_kernel_expression_lhs<ScalarCl<double>>::value));
  EXPECT_FALSE(
      (is_nonscalar_prim_or_rev_kernel_expression<ScalarCl<double>>::value));
  EXPECT_FALSE((is_rev_kernel_expression<ScalarCl<double>>::value));
}

TEST(ScalarCl, value_and_scalar_type) {
  EXPECT_TRUE(
      (std::is_same<stan::scalar_type_t<ScalarCl<double>>, double>::value));
  EXPECT_TRUE((std::is_same<stan::value_type_t<const ScalarCl<double>&>,
                            double>::value));
  EXPECT_FALSE((stan::is_stan_scalar<ScalarCl<double>>::value));
  EXPECT_FALSE((stan::is_autodiff_v<ScalarCl<double>>));
}

TEST(ScalarCl, no_implicit_host_conversions) {
  EXPECT_FALSE((std::is_convertible<ScalarCl<double>, double>::value));
  EXPECT_FALSE((std::is_convertible<double, ScalarCl<double>>::value));
  EXPECT_FALSE((std::is_convertible<int, ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_constructible<ScalarCl<double>, double>::value));
}

TEST(ScalarCl, default_is_zero) {
  ScalarCl<double> x;
  EXPECT_EQ(x.matrix().rows(), 1);
  EXPECT_EQ(x.matrix().cols(), 1);
  EXPECT_EQ(to_host(x), 0.0);
}

TEST(ScalarCl, construct_from_host) {
  ScalarCl<double> x(3.5);
  EXPECT_EQ(to_host(x), 3.5);
  double v = -1.25;
  ScalarCl<double> y(v);
  v = 100.0;
  EXPECT_EQ(to_host(y), -1.25);
}

TEST(ScalarCl, copy_duplicates_device_value) {
  ScalarCl<double> x(2.0);
  ScalarCl<double> y(x);
  EXPECT_NE(x.buffer()(), y.buffer()());
  EXPECT_EQ(to_host(y), 2.0);
  y = ScalarCl<double>(7.0);
  EXPECT_EQ(to_host(x), 2.0);
  EXPECT_EQ(to_host(y), 7.0);

  ScalarCl<double> z;
  z = x;
  EXPECT_NE(x.buffer()(), z.buffer()());
  EXPECT_EQ(to_host(z), 2.0);
}

TEST(ScalarCl, move_transfers_buffer) {
  ScalarCl<double> x(4.0);
  cl_mem x_buf = x.buffer()();
  ScalarCl<double> y(std::move(x));
  EXPECT_EQ(y.buffer()(), x_buf);
  EXPECT_EQ(to_host(y), 4.0);

  ScalarCl<double> z;
  z = std::move(y);
  EXPECT_EQ(z.buffer()(), x_buf);
  EXPECT_EQ(to_host(z), 4.0);
}

#endif
