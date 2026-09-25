#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

using Eigen::MatrixXd;
using stan::math::from_matrix_cl;
using stan::math::matrix_cl;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;

namespace {
MatrixXd test_matrix() {
  MatrixXd m(3, 2);
  m << 1, -2, 3.5, 4, -5, 6;
  return m;
}
}  // namespace

TEST(KernelGeneratorScalarCl, broadcasts_in_matrix_expressions) {
  MatrixXd m = test_matrix();
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(2.5);

  matrix_cl<double> add_cl = m_cl + s;
  EXPECT_MATRIX_EQ(from_matrix_cl(add_cl), (m.array() + 2.5).matrix());
  matrix_cl<double> sub_cl = s - m_cl;
  EXPECT_MATRIX_EQ(from_matrix_cl(sub_cl), (2.5 - m.array()).matrix());
  matrix_cl<double> mul_cl = s * m_cl;
  EXPECT_MATRIX_EQ(from_matrix_cl(mul_cl), (m * 2.5));
  matrix_cl<double> mul2_cl = m_cl * s;
  EXPECT_MATRIX_EQ(from_matrix_cl(mul2_cl), (m * 2.5));
  matrix_cl<double> div_cl = stan::math::elt_divide(m_cl, s);
  EXPECT_MATRIX_EQ(from_matrix_cl(div_cl), (m / 2.5));
  matrix_cl<double> sel_cl = stan::math::select(m_cl > 0.0, s, m_cl);
  MatrixXd sel = (m.array() > 0.0).select(2.5, m);
  EXPECT_MATRIX_EQ(from_matrix_cl(sel_cl), sel);
  MatrixXd m_abs = m.array().abs();
  matrix_cl<double> m_abs_cl(m_abs);
  matrix_cl<double> fn_cl = stan::math::pow(m_abs_cl, s);
  EXPECT_MATRIX_NEAR(from_matrix_cl(fn_cl), m_abs.array().pow(2.5).matrix(),
                     1e-12);
}

TEST(KernelGeneratorScalarCl, same_scalar_used_twice) {
  MatrixXd m = test_matrix();
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(3.0);
  ScalarCl<double> t(-1.0);
  matrix_cl<double> res_cl = m_cl * s + s - t;
  EXPECT_MATRIX_EQ(from_matrix_cl(res_cl),
                   (m.array() * 3.0 + 3.0 + 1.0).matrix());
}

TEST(KernelGeneratorScalarCl, matrix_compound_assign_with_scalar) {
  MatrixXd m = test_matrix();
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(0.5);
  m_cl += s;
  EXPECT_MATRIX_EQ(from_matrix_cl(m_cl), (m.array() + 0.5).matrix());
}

TEST(KernelGeneratorScalarCl, scalar_only_operators_are_eager) {
  ScalarCl<double> a(3.0);
  ScalarCl<double> b(-2.0);
  EXPECT_TRUE((std::is_same<decltype(a + b), ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(a * 2.0), ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(-a), ScalarCl<double>>::value));
  EXPECT_EQ(to_host(a + b), 1.0);
  EXPECT_EQ(to_host(a - b), 5.0);
  EXPECT_EQ(to_host(a * b), -6.0);
  EXPECT_EQ(to_host(a / b), -1.5);
  EXPECT_EQ(to_host(a + 1.0), 4.0);
  EXPECT_EQ(to_host(1.0 - a), -2.0);
  EXPECT_EQ(to_host(2 * a), 6.0);
  EXPECT_EQ(to_host(a / 4.0), 0.75);
  EXPECT_EQ(to_host(-b), 2.0);
}

TEST(KernelGeneratorScalarCl, scalar_only_functions_are_eager) {
  ScalarCl<double> a(0.25);
  ScalarCl<double> b(2.0);
  EXPECT_TRUE((std::is_same<decltype(exp(a)), ScalarCl<double>>::value));
  EXPECT_NEAR(to_host(exp(a)), std::exp(0.25), 1e-15);
  EXPECT_NEAR(to_host(log(a)), std::log(0.25), 1e-15);
  EXPECT_NEAR(to_host(sqrt(a)), 0.5, 1e-15);
  EXPECT_NEAR(to_host(square(b)), 4.0, 1e-15);
  EXPECT_NEAR(to_host(stan::math::opencl::inv_logit(a)),
              1.0 / (1.0 + std::exp(-0.25)), 1e-15);
  EXPECT_NEAR(to_host(pow(a, b)), 0.0625, 1e-15);
  EXPECT_NEAR(to_host(pow(b, 3.0)), 8.0, 1e-15);
  EXPECT_NEAR(to_host(fmax(a, b)), 2.0, 1e-15);
  EXPECT_NEAR(to_host(ldexp(b, 3)), 16.0, 1e-15);
  EXPECT_NEAR(to_host(exp(a) * b + log(b)),
              2.0 * std::exp(0.25) + std::log(2.0), 1e-15);
}

TEST(KernelGeneratorScalarCl, construct_from_expression) {
  MatrixXd one(1, 1);
  one << 4.0;
  matrix_cl<double> one_cl(one);
  ScalarCl<double> s(2.0);
  ScalarCl<double> r(one_cl * s + 1.0);
  EXPECT_EQ(to_host(r), 9.0);
  r = one_cl - s;
  EXPECT_EQ(to_host(r), 2.0);
}

TEST(KernelGeneratorScalarCl, rejects_matrix_sized_expressions) {
  matrix_cl<double> m_cl(test_matrix());
  matrix_cl<double> empty_cl(0, 0);
  ScalarCl<double> s(2.0);
  EXPECT_THROW(ScalarCl<double>(m_cl * s), std::invalid_argument);
  EXPECT_THROW(ScalarCl<double>(m_cl + 1.0), std::invalid_argument);
  EXPECT_THROW(s = m_cl * 2.0, std::invalid_argument);
  EXPECT_THROW(s += m_cl, std::invalid_argument);
  EXPECT_THROW(ScalarCl<double>(empty_cl + 1.0), std::invalid_argument);
  EXPECT_EQ(to_host(s), 2.0);
}

TEST(KernelGeneratorScalarCl, compound_assignment) {
  ScalarCl<double> s(2.0);
  ScalarCl<double> t(3.0);
  s += t;
  EXPECT_EQ(to_host(s), 5.0);
  s -= 1.0;
  EXPECT_EQ(to_host(s), 4.0);
  s *= t;
  EXPECT_EQ(to_host(s), 12.0);
  s /= 4.0;
  EXPECT_EQ(to_host(s), 3.0);
  MatrixXd one(1, 1);
  one << 10.0;
  matrix_cl<double> one_cl(one);
  s += one_cl;
  EXPECT_EQ(to_host(s), 13.0);
  EXPECT_EQ(to_host(t), 3.0);
}

TEST(KernelGeneratorScalarCl, chained_scalar_results_stay_ordered) {
  MatrixXd m = test_matrix();
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(1.0);
  for (int i = 0; i < 20; ++i) {
    s = s * 1.5 + 0.25;
  }
  double expected = 1.0;
  for (int i = 0; i < 20; ++i) {
    expected = expected * 1.5 + 0.25;
  }
  matrix_cl<double> res_cl = m_cl * s;
  EXPECT_MATRIX_NEAR(from_matrix_cl(res_cl), m * expected, 1e-9);
}

TEST(KernelGeneratorScalarCl, scalar_only_subexpressions_keep_dense_view) {
  // Scalar-only subexpressions have dynamic size; they must not shrink the
  // triangular view inferred for the result.
  MatrixXd m = test_matrix();
  matrix_cl<double> m_cl(m);
  ScalarCl<double> s(2.0);
  auto s_op = stan::math::as_operation_cl(s);
  matrix_cl<double> div_cl
      = stan::math::elt_multiply(stan::math::elt_divide(1.0, s_op), m_cl);
  EXPECT_EQ(div_cl.view(), stan::math::matrix_cl_view::Entire);
  EXPECT_MATRIX_EQ(from_matrix_cl(div_cl), m / 2.0);
  matrix_cl<double> log_cl
      = stan::math::elt_multiply(stan::math::log(s_op), m_cl);
  EXPECT_EQ(log_cl.view(), stan::math::matrix_cl_view::Entire);
  EXPECT_MATRIX_NEAR(from_matrix_cl(log_cl), m * std::log(2.0), 1e-14);
  matrix_cl<double> fmax_cl
      = stan::math::elt_multiply(stan::math::fmax(s_op, 1.0), m_cl);
  EXPECT_EQ(fmax_cl.view(), stan::math::matrix_cl_view::Entire);
  EXPECT_MATRIX_EQ(from_matrix_cl(fmax_cl), m * 2.0);
}

namespace {
/**
 * Checks generated kernel source against the stored reference kernel.
 * Defining STAN_TEST_KERNEL_GENERATOR_STORE_REFERENCE_KERNELS rewrites the
 * reference instead.
 */
void expect_reference_kernel(const std::string& kernel_filename,
                             const std::string& kernel_src) {
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);
}
}  // namespace

TEST(KernelGeneratorScalarCl, broadcast_kernel_source) {
  // the device scalar is read from element 0 of its buffer, and a scalar used
  // twice is passed to the kernel once
  matrix_cl<double> m_cl(test_matrix());
  ScalarCl<double> s(2.0);
  auto expr = m_cl * s + s;
  matrix_cl<double> res_cl;
  expect_reference_kernel("scalar_cl_broadcast.cl",
                          expr.get_kernel_source_for_evaluating_into(res_cl));
  res_cl = expr;
  EXPECT_MATRIX_EQ(from_matrix_cl(res_cl),
                   (test_matrix().array() * 2.0 + 2.0).matrix());
}

TEST(KernelGeneratorScalarCl, scalar_result_kernel_source) {
  // a scalar-only expression is evaluated by one thread into a 1x1 buffer
  ScalarCl<double> s(2.0);
  ScalarCl<double> t(3.0);
  auto expr = stan::math::scalar_result_<decltype(
      stan::math::exp(stan::math::as_operation_cl(s))
      + stan::math::as_operation_cl(t))>(
      stan::math::exp(stan::math::as_operation_cl(s))
      + stan::math::as_operation_cl(t));
  ScalarCl<double> res;
  expect_reference_kernel(
      "scalar_cl_scalar_result.cl",
      expr.get_kernel_source_for_evaluating_into(res.matrix()));
  res.matrix() = expr;
  EXPECT_NEAR(to_host(res), std::exp(2.0) + 3.0, 1e-14);
}

TEST(KernelGeneratorScalarCl, compound_assignment_kernel_source) {
  // s += t reads and writes the same buffer in one single-thread kernel
  ScalarCl<double> s(2.0);
  ScalarCl<double> t(3.0);
  auto expr = stan::math::scalar_result_<decltype(
      stan::math::as_operation_cl(s) + stan::math::as_operation_cl(t))>(
      stan::math::as_operation_cl(s) + stan::math::as_operation_cl(t));
  expect_reference_kernel(
      "scalar_cl_compound_assignment.cl",
      expr.get_kernel_source_for_evaluating_into(s.matrix()));
  s += t;
  EXPECT_EQ(to_host(s), 5.0);
}

TEST(KernelGeneratorScalarCl, fused_check_kernel_source) {
  // a check of a device scalar placed first in a multi-result kernel: the
  // number of threads comes from the matrix expression
  MatrixXd m = test_matrix();
  matrix_cl<double> m_cl(m);
  ScalarCl<double> sigma(1.5);
  auto sigma_op = stan::math::as_operation_cl(sigma);
  auto check_sigma
      = stan::math::check_cl("test", "sigma", sigma_op, "positive");
  matrix_cl<double> res_cl;
  auto res = stan::math::results(check_sigma, res_cl);
  auto exprs = stan::math::expressions(sigma_op > 0.0,
                                       stan::math::elt_divide(m_cl, sigma_op));
  expect_reference_kernel("scalar_cl_fused_check.cl",
                          res.get_kernel_source_for_evaluating(exprs));
  res = exprs;
  EXPECT_MATRIX_NEAR(from_matrix_cl(res_cl), m / 1.5, 1e-14);
}

#endif
