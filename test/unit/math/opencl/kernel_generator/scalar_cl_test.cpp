#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
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

#endif
