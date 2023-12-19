#ifndef TEST_UNIT_MATH_TEST_AD_HPP
#define TEST_UNIT_MATH_TEST_AD_HPP

#include <stan/math/mix.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/ad_tolerances.hpp>
#include <test/unit/math/is_finite.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/test_ad_matvar.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using d_t = double;
using v_t = stan::math::var;
using fd_t = stan::math::fvar<d_t>;
using ffd_t = stan::math::fvar<fd_t>;
using fv_t = stan::math::fvar<stan::math::var>;
using ffv_t = stan::math::fvar<fv_t>;

namespace stan {
namespace test {
namespace internal {

/**
 * Evaluates nested matrix template expressions, which is a no-op for
 * arithmetic arguments.
 *
 * @tparam T arithmetic type
 * @param[in] x value
 * @return value
 */
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
auto eval(T x) {
  return x;
}

/**
 * Evaluates nested matrix template expressions, which is a no-op for
 * complex arguments.
 *
 * @tparam T complex value type
 * @param x[in] value
 * @return value
 */
template <typename T>
auto eval(const std::complex<T>& x) {
  return x;
}

/**
 * Evaluates all nested matrix expression templates, which is a no-op for
 * reverse-mode autodiff variables.
 *
 * @param[in] x value
 * @return value
 */
auto eval(const stan::math::var& x) { return x; }

/**
 * Evaluates all matrix expression templates, which is a no-op for
 * forward-mode autodiff variables.
 *
 * @tparam T value type of fvar
 * @param[in] x value
 * @return value
 */
template <typename T>
auto eval(const stan::math::fvar<T>& x) {
  return x;
}

/**
 * Evaluates all nested matrix expression templates, which evaluates
 * the specified derived matrix.
 *
 * @tparam Derived derived type of the expression
 * @param x expression
 * @return evaluated expression
 */
template <typename Derived>
auto eval(const Eigen::EigenBase<Derived>& x) {
  return x.derived().eval();
}

/**
 * Evaluates all nested matrix expression templates elementwise.
 *
 * @tparam T type of elements
 * @param[in] x vector of expressions
 * @return vector of evaluated expressions
 */
template <typename T>
auto eval(const std::vector<T>& x) {
  using T_res = decltype(eval(std::declval<T>()));
  std::vector<T_res> res;
  for (auto& i : x) {
    res.push_back(eval(i));
  }
  return res;
}

/**
 * Tests that the specified function applied to the specified argument
 * yields the values and gradients consistent with finite differences
 * applied to the `double` values.  Gradients are calculated using the
 * functional `stan::math::gradient<F>`, which uses autodiff variables
 * of type `stan::math::var.`
 *
 * <p>All tests are done using relative values and accounting for
 * infinite, not-a-number, and zero values, as defined in the
 * documentation for `expect_near`.
 *
 * @tparam F type of functor
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_gradient(const ad_tolerances& tols, const F& f,
                   const Eigen::VectorXd& x, double fx,
                   bool test_derivs = true) {
  Eigen::VectorXd grad_ad;
  double fx_ad = fx;
  stan::math::gradient<F>(f, x, fx_ad, grad_ad);
  expect_near_rel("gradient() val", fx, fx_ad, tols.gradient_val_);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  stan::math::finite_diff_gradient_auto(f, x, fx_fd, grad_fd);
  expect_near_rel("gradient() grad for finite diff vs auto diff", grad_fd,
                  grad_ad, tols.gradient_grad_);
}

#ifndef STAN_MATH_TESTS_REV_ONLY
/**
 * Tests that the specified function applied to the specified argument
 * yields the values and gradients consistent with finite differences
 * applied to the `double` values.  Gradients are calculated using the
 * functional `stan::math::gradient<double, F>`, which uses autodiff
 * variables of type `stan::math::fvar<double>.`
 *
 * <p>The functor to test must define `T operator()(const
 * Matrix<T, -1, 1>& x) const;` for relevant `double` and autodiff
 * types.
 *
 *
 * <p>All tests are done using relative values and accounting for
 * infinite, not-a-number, and zero values, as defined in the
 * documentation for `expect_near`.
 *
 * @tparam F type of functor
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_gradient_fvar(const ad_tolerances& tols, const F& f,
                        const Eigen::VectorXd& x, double fx,
                        bool test_derivs = true) {
  Eigen::VectorXd grad_ad;
  double fx_ad = fx;
  stan::math::gradient<double, F>(f, x, fx_ad, grad_ad);
  expect_near_rel("gradient_fvar() val", fx, fx_ad, tols.gradient_fvar_val_);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  stan::math::finite_diff_gradient_auto(f, x, fx_fd, grad_fd);
  expect_near_rel("gradient_fvar() grad", grad_fd, grad_ad,
                  tols.gradient_fvar_grad_);
}

/**
 * Tests that the specified function applied to the specified argument
 * yields values, gradients, and Hessian consistent with finite
 * differences applied to the `double` values.  The Hessian and
 * gradients are calculated using the functional
 * `stan::math::hessian<double, F>`, which uses autodiff variables of
 * type `stan::math::fvar<fvar<double>>.`
 *
 * <p>The functor to test must define `T operator()(const
 * Matrix<T, -1, 1>& x) const;` for relevant `double` and autodiff
 * types.
 *
 * <p>All tests are done using relative values and accounting for
 * infinite, not-a-number, and zero values, as defined in the
 * documentation for `expect_near`.
 *
 * @tparam F type of functor
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_hessian_fvar(const ad_tolerances& tols, const F& f,
                       const Eigen::VectorXd& x, double fx,
                       bool test_derivs = true) {
  double fx_ad;
  Eigen::VectorXd grad_ad;
  Eigen::MatrixXd H_ad;
  stan::math::hessian<double, F>(f, x, fx_ad, grad_ad, H_ad);
  expect_near_rel("hessian_fvar() val", fx, fx_ad, tols.hessian_fvar_val_);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  stan::math::internal::finite_diff_hessian_auto(f, x, fx_fd, grad_fd, H_fd);
  expect_near_rel("hessian_fvar() grad", grad_fd, grad_ad,
                  tols.hessian_fvar_grad_);
  expect_near_rel("hessian_fvar() Hessian", H_fd, H_ad,
                  tols.hessian_fvar_hessian_);
}

/**
 * Tests that the specified function applied to the specified argument
 * yields values, gradients, and Hessian consistent with finite
 * differences applied to the `double` values.  The Hessian and
 * gradients are calculated using the functional
 * `stan::math::hessian<F>`, which uses autodiff variables of
 * type `stan::math::fvar<var>.`
 *
 * <p>The functor to test must define `T operator()(const
 * Matrix<T, -1, 1>& x) const;` for relevant `double` and autodiff
 * types.
 *
 * <p>All tests are done using relative values and accounting for
 * infinite, not-a-number, and zero values, as defined in the
 * documentation for `expect_near`.
 *
 * @tparam F type of functor
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_hessian(const ad_tolerances& tols, const F& f,
                  const Eigen::VectorXd& x, double fx,
                  bool test_derivs = true) {
  double fx_ad;
  Eigen::VectorXd grad_ad;
  Eigen::MatrixXd H_ad;
  stan::math::hessian<F>(f, x, fx_ad, grad_ad, H_ad);
  expect_near_rel("hessian val", fx, fx_ad, tols.hessian_val_);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  stan::math::internal::finite_diff_hessian_auto(f, x, fx_fd, grad_fd, H_fd);
  expect_near_rel("hessian() grad", grad_fd, grad_ad, tols.hessian_grad_);
  expect_near_rel("hessian() Hessian", H_fd, H_ad, tols.hessian_hessian_);
}

/**
 * Tests that the specified function applied to the specified argument
 * yields values, Hessian, and gradient of Hessian consistent with finite
 * differences applied to the `double` values.  The Hessian and
 * gradient of Hessian are calculated using the functional
 * `stan::math::grad_hessian<F>`, which uses autodiff variables of
 * type `stan::math::fvar<fvar<var>>.`
 *
 * <p>The functor to test must define `T operator()(const
 * Matrix<T, -1, 1>& x) const;` for relevant `double` and autodiff
 * types.
 *
 * <p>All tests are done using relative values and accounting for
 * infinite, not-a-number, and zero values, as defined in the
 * documentation for `expect_near`.
 *
 * @tparam F type of functor
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_grad_hessian(const ad_tolerances& tols, const F& f,
                       const Eigen::VectorXd& x, double fx,
                       bool test_derivs = true) {
  double fx_ad;
  Eigen::MatrixXd H_ad;
  std::vector<Eigen::MatrixXd> grad_H_ad;
  stan::math::grad_hessian(f, x, fx_ad, H_ad, grad_H_ad);
  expect_near_rel("grad_hessian() val", fx, fx_ad, tols.grad_hessian_val_);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::MatrixXd H_fd;
  std::vector<Eigen::MatrixXd> grad_H_fd;
  stan::math::finite_diff_grad_hessian_auto(f, x, fx_fd, H_fd, grad_H_fd);
  expect_near_rel("grad_hessian() Hessian", H_fd, H_ad,
                  tols.grad_hessian_hessian_);
  EXPECT_EQ(x.size(), grad_H_fd.size());
  for (size_t i = 0; i < grad_H_fd.size(); ++i)
    expect_near_rel("grad_hessian() grad Hessian", grad_H_fd[i], grad_H_ad[i],
                    tols.grad_hessian_grad_hessian_);
}
#endif

/**
 * For the specified functor and argument, test that automatic
 * differentiation provides the same value as the double-based version and
 * the same derivatives as finite differences over the double
 * version.  It also tests that an autodiff version throws an
 * exception if and only if the double-based version throws an
 * exception.  It also ensures consistency of infinite, not-a-number,
 * and zero results as defined by the function `expect_near`.
 *
 * <p>The functor to test must define `T operator()(const
 * Matrix<T, -1, 1>& x) const;` for relevant `double` and autodiff types.
 *
 * @tparam G type of polymorphic functor
 * @param tols tolerances for test
 * @param g polymorphic functor from vectors to scalars
 * @param x argument to test
 */
template <typename G>
void expect_ad_derivatives(const ad_tolerances& tols, const G& g,
                           const Eigen::VectorXd& x) {
  double gx = g(x);
  test_gradient(tols, g, x, gx);
#ifndef STAN_MATH_TESTS_REV_ONLY
  test_gradient_fvar(tols, g, x, gx);
  test_hessian(tols, g, x, gx);
  test_hessian_fvar(tols, g, x, gx);
  test_grad_hessian(tols, g, x, gx);
#endif
}

/**
 * Test that the specified functor applied to the specified value
 * throws an exception, reporting the name of the expected exception
 * type as provided in the case of failure and raising a google test
 * error.
 *
 * @tparam T type of scalars to test
 * @tparam F type of functor to test
 * @param f functor to test
 * @param x value to test
 * @param name_of_T name of type of exception expected
 */
template <typename T, typename F>
void expect_throw(const F& f, const Eigen::VectorXd& x,
                  const std::string& name_of_T) {
  Eigen::Matrix<T, -1, 1> x_t(x.rows());
  for (int i = 0; i < x.rows(); ++i)
    x_t(i) = x(i);
  try {
    f(x_t);
    FAIL() << "double throws, expect type " << name_of_T
           << " version to throw for x = " << x;
  } catch (...) {
    SUCCEED();
  }
}

/**
 * Succeeds if the specified function applied to the specified
 * argument throws an exception at every level of autodiff.
 *
 * @tparam F type of function
 * @param f function to test
 * @param x argument to test
 */
template <typename F>
void expect_all_throw(const F& f, const Eigen::VectorXd& x) {
  using stan::math::fvar;
  using stan::math::var;
  expect_throw<double>(f, x, "double");
  expect_throw<var>(f, x, "var");
#ifndef STAN_MATH_TESTS_REV_ONLY
  expect_throw<fvar<double>>(f, x, "fvar<double>");
  expect_throw<fvar<fvar<double>>>(f, x, "fvar<fvar<double>>");
  expect_throw<fvar<var>>(f, x, "fvar<var>");
  expect_throw<fvar<fvar<var>>>(f, x, "fvar<fvar<var>>");
#endif
}

/**
 * Succeeds if the specified function applied to the specified
 * argument throws an exception at every level of autodiff.
 *
 * @tparam F type of function
 * @param f function to evaluate
 * @param x argument to evaluate
 */
template <typename F>
void expect_all_throw(const F& f, double x1) {
  using stan::math::serialize_return;
  auto h = [&](auto v) { return serialize_return(eval(f(v(0)))); };
  Eigen::VectorXd x(1);
  x << x1;
  expect_all_throw(h, x);
}

/**
 * Succeeds if the specified function applied to the specified
 * argument throws an exception at every level of autodiff.
 *
 * @tparam F type of function
 * @param f function to evaluate
 * @param x1 first argument
 * @param x2 second argument
 */
template <typename F>
void expect_all_throw(const F& f, double x1, double x2) {
  using stan::math::serialize_return;
  auto h = [&](auto v) { return serialize_return(eval(f(v(0), v(1)))); };
  Eigen::VectorXd x(2);
  x << x1, x2;
  expect_all_throw(h, x);
}

/**
 * Succeeds if the specified function applied to the specified
 * argument throws an exception at every level of autodiff.
 *
 * @tparam F type of function
 * @param f function to evaluate
 * @param x1 first argument
 * @param x2 second argument
 * @param x3 third argument
 */
template <typename F>
void expect_all_throw(const F& f, double x1, double x2, double x3) {
  using stan::math::serialize_return;
  auto h = [&](auto v) { return serialize_return(eval(f(v(0), v(1), v(2)))); };
  Eigen::VectorXd x(3);
  x << x1, x2, x3;
  expect_all_throw(h, x);
}

/**
 * For the specified functor, serialized form of the functor,
 * serialized argument and raw argument sequence, test that automatic
 * differentiation at all levels provides the same answer as the
 * double-based version with finite differences.
 *
 * @tparam F type of original functor applying to sequence of
 * arguments
 * @tparam H type of serialized functor applying to Eigen vector and
 * returning a single component of the value
 * @tparam Ts type pack for arguments to original functor with double
 * scalar types
 * @param tols tolerances for test
 * @param f functor to evaluate
 * @param g serialized functor taking an Eigen vector and returning a
 * serialized container of the original output
 * that returns that component of the original output
 * original output
 * @param x serialized input
 * @param xs sequence of arguments with double-based scalars
 */
template <typename F, typename G, typename... Ts>
void expect_ad_helper(const ad_tolerances& tols, const F& f, const G& g,
                      const Eigen::VectorXd& x, Ts... xs) {
  using stan::math::serialize;
  auto h
      = [&](const int i) { return [&g, i](const auto& v) { return g(v)[i]; }; };
  size_t result_size = 0;
  try {
    auto y1 = eval(f(xs...));  // original types, including int
    auto y2 = eval(g(x));      // all int cast to double
    auto y1_serial = serialize<double>(y1);
    expect_near_rel("expect_ad_helper", y1_serial, y2, 1e-10);
    result_size = y1_serial.size();
  } catch (...) {
    internal::expect_all_throw(h(0), x);
    return;
  }
  for (size_t i = 0; i < result_size; ++i) {
    expect_ad_derivatives(tols, h(i), x);
  }
}

/**
 * Test that the specified unary functor and arguments produce for
 * every autodiff type the same value as the specified argument and
 * the same derivatives as finite differences.
 *
 * @tparam F type of functor to test
 * @tparam T type of first argument with double-based scalar
 * @param tols tolerances for test
 * @param f functor to test
 * @param x argument to test
 */
template <typename F, typename T>
void expect_ad_v(const ad_tolerances& tols, const F& f, const T& x) {
  using stan::math::serialize_args;
  using stan::math::serialize_return;
  using stan::math::to_deserializer;
  auto g = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto xds = ds.read(x);
    return serialize_return(eval(f(xds)));
  };
  internal::expect_ad_helper(tols, f, g, serialize_args(x), x);
}

/**
 * Test that the specified unary functor and arguments produce for
 * every autodiff type the same value as the specified argument and
 * the same derivatives as finite differences.
 *
 * This is an overload of this function for integer arguments to make
 * sure that if the integer version throws, so do all the autodiff
 * versions and to test that autodiff works at the integer specified
 * cast to a double.
 *
 * @tparam F type of functor to test
 * @tparam T type of first argument with double-based scalar
 * @param tols tolerances for test
 * @param f functor to test
 * @param x argument to test
 */
template <typename F>
void expect_ad_v(const ad_tolerances& tols, const F& f, int x) {
  double x_dbl = static_cast<double>(x);

  // if f throws on int, must throw everywhere with double
  try {
    f(x);
  } catch (...) {
    expect_all_throw(f, x_dbl);
    return;
  }

  // values must be same with int and double
  expect_near_rel("expect_ad_v: int must produce same result as double cast",
                  f(x_dbl), f(x));

  // autodiff should work at double value
  expect_ad_v(tols, f, x_dbl);
}

/**
 * Test that the specified binary functor and arguments produce for
 * every autodiff type the same value as the double-based version and
 * the same derivatives as finite differences when both arguments are
 * autodiff variables.
 *
 * @tparam F type of functor to test
 * @tparam T1 type of first argument with double-based scalar
 * @tparam T2 type of second argument with double-based scalar
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 first argument
 * @param x2 second argument
 */
template <typename F, typename T1, typename T2>
void expect_ad_vv(const ad_tolerances& tols, const F& f, const T1& x1,
                  const T2& x2) {
  using stan::math::serialize_args;
  using stan::math::serialize_return;
  using stan::math::to_deserializer;
  // d.x1
  auto g1 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    return serialize_return(eval(f(x1ds, x2)));
  };
  internal::expect_ad_helper(tols, f, g1, serialize_args(x1), x1, x2);

  // d.x2
  auto g2 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x2ds = ds.read(x2);
    return serialize_return(eval(f(x1, x2ds)));
  };
  internal::expect_ad_helper(tols, f, g2, serialize_args(x2), x1, x2);

  // d.x1, d.x2
  auto g12 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    auto x2ds = ds.read(x2);
    return serialize_return(eval(f(x1ds, x2ds)));
  };
  internal::expect_ad_helper(tols, f, g12, serialize_args(x1, x2), x1, x2);
}

template <typename F, typename T2>
void expect_ad_vv(const ad_tolerances& tols, const F& f, int x1, const T2& x2) {
  try {
    f(x1, x2);
  } catch (...) {
    expect_all_throw(f, x1, x2);
    return;
  }

  double x1_dbl = static_cast<double>(x1);

  // expect same result with int or and cast to double
  expect_near_rel("expect_ad_vv(int, T2)", f(x1, x2), f(x1_dbl, x2));

  // expect autodiff to work at double value
  expect_ad_vv(tols, f, x1_dbl, x2);

  // expect autodiff to work when binding int; includes expect-all-throw test
  auto g = [&](const auto& u) { return f(x1, u); };
  expect_ad_v(tols, g, x2);
}

template <typename F, typename T1>
void expect_ad_vv(const ad_tolerances& tols, const F& f, const T1& x1, int x2) {
  try {
    f(x1, x2);
  } catch (...) {
    expect_all_throw(f, x1, x2);
    return;
  }

  double x2_dbl = static_cast<double>(x2);

  // expect same result with int or and cast to double
  expect_near_rel("expect_ad_vv(T1, int)", f(x1, x2), f(x1, x2_dbl));

  // expect autodiff to work at double value
  expect_ad_vv(tols, f, x1, x2_dbl);

  // expect autodiff to work when binding int; includes expect-all-throw test
  auto g = [&](const auto& u) { return f(u, x2); };
  expect_ad_v(tols, g, x1);
}

template <typename F>
void expect_ad_vv(const ad_tolerances& tols, const F& f, int x1, int x2) {
  // this one needs throw test because it's not handled by recursion
  try {
    f(x1, x2);
  } catch (...) {
    expect_all_throw(f, x1, x2);
    return;
  }

  double x1_dbl = static_cast<double>(x1);
  double x2_dbl = static_cast<double>(x2);

  // expect same result with int or and cast to double
  expect_near_rel("expect_ad_vv(int, int)", f(x1, x2), f(x1_dbl, x2_dbl));

  // expect autodiff to work at double values
  // these take care of x1_dbl, x2_dbl case by delegation
  // they also take care of binding int tests
  expect_ad_vv(tols, f, x1, x2_dbl);
  expect_ad_vv(tols, f, x1_dbl, x2);
}

/**
 * Test that the specified ternary functor and arguments produce for
 * every autodiff type the same value as the double-based version and
 * the same derivatives as finite differences when both arguments are
 * autodiff variables.
 *
 * @tparam F type of functor to test
 * @tparam T1 type of first argument with double-based scalar
 * @tparam T2 type of second argument with double-based scalar
 * @tparam T3 type of third argument with double-based scalar
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 first argument
 * @param x2 second argument
 * @param x3 third argument
 */
template <typename F, typename T1, typename T2, typename T3>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, const T1& x1,
                   const T2& x2, const T3& x3) {
  using stan::math::serialize_args;
  using stan::math::serialize_return;
  using stan::math::to_deserializer;
  // d.x1
  auto g1 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    return serialize_return(eval(f(x1ds, x2, x3)));
  };
  internal::expect_ad_helper(tols, f, g1, serialize_args(x1), x1, x2, x3);

  // d.x2
  auto g2 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x2ds = ds.read(x2);
    return serialize_return(eval(f(x1, x2ds, x3)));
  };
  internal::expect_ad_helper(tols, f, g2, serialize_args(x2), x1, x2, x3);

  // d.x3
  auto g3 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x3ds = ds.read(x3);
    return serialize_return(eval(f(x1, x2, x3ds)));
  };
  internal::expect_ad_helper(tols, f, g3, serialize_args(x3), x1, x2, x3);

  // d.x1 d.x2
  auto g12 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    auto x2ds = ds.read(x2);
    return serialize_return(eval(f(x1ds, x2ds, x3)));
  };
  internal::expect_ad_helper(tols, f, g12, serialize_args(x1, x2), x1, x2, x3);

  // d.x1 d.x3
  auto g13 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    auto x3ds = ds.read(x3);
    return serialize_return(eval(f(x1ds, x2, x3ds)));
  };
  internal::expect_ad_helper(tols, f, g13, serialize_args(x1, x3), x1, x2, x3);

  // d.x2 d.x3
  auto g23 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x2ds = ds.read(x2);
    auto x3ds = ds.read(x3);
    return serialize_return(eval(f(x1, x2ds, x3ds)));
  };
  internal::expect_ad_helper(tols, f, g23, serialize_args(x2, x3), x1, x2, x3);

  // d.x1 d.x2 d.x3
  auto g123 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    auto x2ds = ds.read(x2);
    auto x3ds = ds.read(x3);
    return serialize_return(eval(f(x1ds, x2ds, x3ds)));
  };
  internal::expect_ad_helper(tols, f, g123, serialize_args(x1, x2, x3), x1, x2,
                             x3);
}

template <typename F, typename T3>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, int x1, int x2,
                   const T3& x3) {
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x1_dbl = static_cast<double>(x1);
  double x2_dbl = static_cast<double>(x2);

  // test all promotion patterns;  includes x1_dbl & x2_dbl recursively
  expect_ad_vvv(tols, f, x1_dbl, x2, x3);
  expect_ad_vvv(tols, f, x1, x2_dbl, x3);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1_dbl, x2_dbl, x3));

  // bind ints and test autodiff
  auto g23 = [=](const auto& u2, const auto& u3) { return f(x1, u2, u3); };
  expect_ad_vv(tols, g23, x2, x3);

  auto g13 = [=](const auto& u1, const auto& u3) { return f(u1, x2, u3); };
  expect_ad_vv(tols, g13, x1, x3);
}

template <typename F, typename T2, typename T3>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, int x1, const T2& x2,
                   const T3& x3) {
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x1_dbl = static_cast<double>(x1);

  // test all promotion patterns
  expect_ad_vvv(tols, f, x1_dbl, x2, x3);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1_dbl, x2, x3));

  // bind ints and test autodiff
  auto g23 = [=](const auto& u2, const auto& u3) { return f(x1, u2, u3); };
  expect_ad_vv(tols, g23, x2, x3);
}

template <typename F, typename T1, typename T3>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, const T1& x1, int x2,
                   const T3& x3) {
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x2_dbl = static_cast<double>(x2);

  // test promotion
  expect_ad_vvv(tols, f, x1, x2_dbl, x3);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1, x2_dbl, x3));

  // bind ints and test autodiff
  auto g13 = [=](const auto& u1, const auto& u3) { return f(u1, x2, u3); };
  expect_ad_vv(tols, g13, x1, x3);
}

template <typename F, typename T1, typename T2>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, const T1& x1,
                   const T2& x2, int x3) {
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x3_dbl = static_cast<double>(x3);

  // test promotion
  expect_ad_vvv(tols, f, x1, x2, x3_dbl);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1, x2, x3_dbl));

  // bind ints and test autodiff
  auto g12 = [=](const auto& u1, const auto& u2) { return f(u1, u2, x3); };
  expect_ad_vv(tols, g12, x1, x2);
}

template <typename F, typename T2>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, int x1, const T2& x2,
                   int x3) {
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x1_dbl = static_cast<double>(x1);
  double x3_dbl = static_cast<double>(x3);

  // test promotion recursively
  expect_ad_vvv(tols, f, x1_dbl, x2, x3);
  expect_ad_vvv(tols, f, x1, x2, x3_dbl);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1_dbl, x2, x3_dbl));

  // bind ints and test autodiff
  auto g23 = [=](const auto& u2, const auto& u3) { return f(x1, u2, u3); };
  expect_ad_vv(tols, g23, x2, x3);

  auto g12 = [=](const auto& u1, const auto& u2) { return f(u1, u2, x3); };
  expect_ad_vv(tols, g12, x1, x2);
}

template <typename F, typename T1>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, const T1& x1, int x2,
                   int x3) {
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x2_dbl = static_cast<double>(x2);
  double x3_dbl = static_cast<double>(x3);

  // test promotion recursively
  expect_ad_vvv(tols, f, x1, x2_dbl, x3);
  expect_ad_vvv(tols, f, x1, x2, x3_dbl);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1, x2_dbl, x3_dbl));

  // bind ints and test autodiff
  auto g13 = [=](const auto& u1, const auto& u3) { return f(u1, x2, u3); };
  expect_ad_vv(tols, g13, x1, x3);

  auto g12 = [=](const auto& u1, const auto& u2) { return f(u1, u2, x3); };
  expect_ad_vv(tols, g12, x1, x2);
}

template <typename F>
void expect_ad_vvv(const ad_tolerances& tols, const F& f, int x1, int x2,
                   int x3) {
  // test exception behavior; other exception cases tested recursively
  try {
    f(x1, x2, x3);
  } catch (...) {
    expect_all_throw(f, x1, x2, x3);
    return;
  }

  double x1_dbl = static_cast<double>(x1);
  double x2_dbl = static_cast<double>(x2);
  double x3_dbl = static_cast<double>(x3);

  // test value
  expect_near_rel("expect_ad_vvv(int, int, T3)", f(x1, x2, x3),
                  f(x1_dbl, x2_dbl, x3_dbl));

  // test all promotion patterns;  includes all combos recursively
  expect_ad_vvv(tols, f, x1_dbl, x2, x3);
  expect_ad_vvv(tols, f, x1, x2_dbl, x3);
  expect_ad_vvv(tols, f, x1, x2, x3_dbl);

  // bind ints and test recursively
  auto g12 = [=](const auto& u1, const auto& u2) { return f(u1, u2, x3); };
  expect_ad_vv(tols, g12, x1, x2);

  auto g13 = [=](const auto& u1, const auto& u3) { return f(u1, x2, u3); };
  expect_ad_vv(tols, g13, x1, x3);

  auto g23 = [=](const auto& u2, const auto& u3) { return f(x1, u2, u3); };
  expect_ad_vv(tols, g23, x2, x3);
}

/**
 * Return a sequence of common non-zero arguments.  This includes
 * positive, negative, positive infinite, negative infinity, and
 * not-a-number values, but does not include zero.
 *
 * @return non-zero arguments
 */
const std::vector<double>& common_nonzero_args() {
  static const std::vector<double> common_nz_args{
      -1.3,
      0.49,
      0.99,
      1.01,
      stan::math::positive_infinity(),
      stan::math::negative_infinity(),
      stan::math::not_a_number()};
  return common_nz_args;
}

/**
 * Return the sequence of common scalar arguments to test.  These
 * include the values returned by `common_nonzero_args()` and zero.
 *
 * @return sequence of common scalar arguments to test
 */
std::vector<double> common_args() {
  auto result = common_nonzero_args();
  result.push_back(0);
  return result;
}

std::vector<int> common_nonzero_int_args() {
  static const std::vector<int> args{-1, 1};
  return args;
}

std::vector<int> common_int_args() {
  std::vector<int> args = common_nonzero_int_args();
  args.push_back(0);
  return args;
}

/**
 * Test that the specified comparison function produces the same result when
 * applied to the specified double or integer values as it does when
 * applied to their promotion to autodiff variables.  The return type
 * for comparisons is boolean, so there are no derivatives to test.
 *
 * <p>Arguments are tested at autodiff levels: reverse (`var`),
 * forward (`fvar<double>`), forward within forward
 * (`fvar<fvar<double>>`), reverse within forward (`fvar<var>`), and
 * reverse within forward within forward (`fvar<fvar<var>>`).  The
 * specified functor must be overloaded to handle all of these
 * possibilities in the first and/or second argument positions.
 *
 * <p>The values being tested (of types `T1` and `T2`) must be `int`
 * or `double`.
 *
 * @tparam F type of polymorphic functor being tested
 * @tparam T1 type of first argument being tested
 * @tparam T2 type of second argument being tested
 * @param f functor being tested
 * @param x1 first value being tested
 * @param x2 second value being tested
 */
template <typename F, typename T1, typename T2>
void expect_comparison(const F& f, const T1& x1, const T2& x2) {
  using stan::math::fvar;
  using stan::math::var;
  typedef var v;
  typedef fvar<double> fd;
  typedef fvar<fvar<double>> ffd;
  typedef fvar<var> fv;
  typedef fvar<fvar<var>> ffv;

  // vv
  EXPECT_EQ(f(x1, x2), f(v(x1), v(x2)));
  EXPECT_EQ(f(x1, x2), f(fd(x1), fd(x2)));
  EXPECT_EQ(f(x1, x2), f(ffd(x1), ffd(x2)));
  EXPECT_EQ(f(x1, x2), f(fv(x1), fv(x2)));
  EXPECT_EQ(f(x1, x2), f(ffv(x1), ffv(x2)));

  // vd
  EXPECT_EQ(f(x1, x2), f(v(x1), x2));
  EXPECT_EQ(f(x1, x2), f(fd(x1), x2));
  EXPECT_EQ(f(x1, x2), f(ffd(x1), x2));
  EXPECT_EQ(f(x1, x2), f(fv(x1), x2));
  EXPECT_EQ(f(x1, x2), f(ffv(x1), x2));

  // dv
  EXPECT_EQ(f(x1, x2), f(x1, v(x2)));
  EXPECT_EQ(f(x1, x2), f(x1, fd(x2)));
  EXPECT_EQ(f(x1, x2), f(x1, ffd(x2)));
  EXPECT_EQ(f(x1, x2), f(x1, fv(x2)));
  EXPECT_EQ(f(x1, x2), f(x1, ffv(x2)));
}

}  // namespace internal

template <typename F>
void expect_all_throw(const F& f, double x) {
  internal::expect_all_throw(f, x);
}

template <typename F>
void expect_all_throw(const F& f, double x1, double x2) {
  internal::expect_all_throw(f, x1, x2);
}

template <typename F>
void expect_all_throw(const F& f, double x1, double x2, double x3) {
  internal::expect_all_throw(f, x1, x2, x3);
}

/**
 * Test that the specified function has the same value when applied to
 * type `T` as it does promoting `T` to all of the autodiff types
 * (`var`, `fvar<double>`, `fvar<fvar<double>>`, `fvar<var>`,
 * `fvar<fvar<var>>`).
 *
 * @tparam F type of function to test
 * @tparam T type of scalar argument
 * @param f function to test
 * @param x value to test
 */
template <typename F, typename T>
void expect_value(const F& f, const T& x) {
  using stan::math::fvar;
  using stan::math::var;
  typedef var v;
  typedef fvar<double> fd;
  typedef fvar<fvar<double>> ffd;
  typedef fvar<var> fv;
  typedef fvar<fvar<var>> ffv;

  double fx = f(x);
  EXPECT_FLOAT_EQ(fx, f(v(x)).val());
  EXPECT_FLOAT_EQ(fx, f(fd(x)).val());
  EXPECT_FLOAT_EQ(fx, f(ffd(x)).val().val());
  EXPECT_FLOAT_EQ(fx, f(fv(x)).val().val());
  EXPECT_FLOAT_EQ(fx, f(ffv(x)).val().val().val());
}

/**
 * Test that the specified function has the same value when applied to
 * type `T1` and `T2` to all of the autodiff types
 * (`var`, `fvar<double>`, `fvar<fvar<double>>`, `fvar<var>`,
 * `fvar<fvar<var>>`).
 *
 * @tparam F type of function to test
 * @tparam T1 type of first scalar argument
 * @tparam T2 type of second scalar argument
 * @param f function to test
 * @param x1 first argument to test
 * @param x1 second argument to test
 */
template <typename F, typename T1, typename T2>
void expect_value(const F& f, const T1& x1, const T2& x2) {
  using stan::math::fvar;
  using stan::math::var;
  typedef var v;
  typedef fvar<double> fd;
  typedef fvar<fvar<double>> ffd;
  typedef fvar<var> fv;
  typedef fvar<fvar<var>> ffv;
  double fx = f(x1, x2);

  // vv
  EXPECT_FLOAT_EQ(fx, f(v(x1), v(x2)).val());
  EXPECT_FLOAT_EQ(fx, f(fd(x1), fd(x2)).val());
  EXPECT_FLOAT_EQ(fx, f(ffd(x1), ffd(x2)).val().val());
  EXPECT_FLOAT_EQ(fx, f(fv(x1), fv(x2)).val().val());
  EXPECT_FLOAT_EQ(fx, f(ffv(x1), ffv(x2)).val().val().val());

  // vd
  EXPECT_FLOAT_EQ(fx, f(v(x1), x2).val());
  EXPECT_FLOAT_EQ(fx, f(fd(x1), x2).val());
  EXPECT_FLOAT_EQ(fx, f(ffd(x1), x2).val().val());
  EXPECT_FLOAT_EQ(fx, f(fv(x1), x2).val().val());
  EXPECT_FLOAT_EQ(fx, f(ffv(x1), x2).val().val().val());

  // dv
  EXPECT_FLOAT_EQ(fx, f(x1, v(x2)).val());
  EXPECT_FLOAT_EQ(fx, f(x1, fd(x2)).val());
  EXPECT_FLOAT_EQ(fx, f(x1, ffd(x2)).val().val());
  EXPECT_FLOAT_EQ(fx, f(x1, fv(x2)).val().val());
  EXPECT_FLOAT_EQ(fx, f(x1, ffv(x2)).val().val().val());
}

/**
 * Test that the specified polymorphic unary functor produces autodiff
 * results consistent with values determined by double or integer
 * inputs and 1st-, 2nd-, and 3rd-order derivatives consistent with
 * finite differences of double inputs at the specified tolerances.
 *
 * @tparam F type of functor to test
 * @tparam T type of argument
 * @param tols tolerances for test
 * @param f function to test
 * @param x argument to test
 */
template <typename F, typename T>
void expect_ad(const ad_tolerances& tols, const F& f, const T& x) {
  internal::expect_ad_v(tols, f, x);
}

/**
 * Test that the specified polymorphic unary functor produces autodiff
 * results consistent with values determined by double or integer
 * inputs and 1st-, 2nd-, and 3rd-order derivatives consistent with
 * finite differences of double inputs at default tolerances.
 *
 * @tparam F type of functor to test
 * @tparam T type of argument
 * @param f function to test
 * @param x argument to test
 */
template <typename F, typename T>
void expect_ad(const F& f, const T& x) {
  ad_tolerances tols;
  expect_ad(tols, f, x);
}

/**
 * Test that the specified binary function produces autodiff values
 * and 1st-, 2nd-, and 3rd-order derivatives consistent with primitive
 * int and double inputs and finite differences at the specified
 * tolerances.
 *
 * @tparam F type of binary polymorphic functor to test
 * @tparam T1 type of double- or int-based first argument
 * @tparam T2 type of double- or int-based second argument
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 first argument to test
 * @param x2 second argument to test
 */
template <typename F, typename T1, typename T2>
void expect_ad(const ad_tolerances& tols, const F& f, const T1& x1,
               const T2& x2) {
  internal::expect_ad_vv(tols, f, x1, x2);
}

/**
 * Test that the specified binary function produces autodiff values
 * and derivatives consistent with primitive int and double inputs and
 * finite differences at default tolerances.
 *
 * @tparam F type of binary polymorphic functor to test
 * @tparam T1 type of double- or int-based first argument
 * @tparam T2 type of double- or int-based second argument
 * @param f functor to test
 * @param x1 first argument to test
 * @param x2 second argument to test
 */
template <typename F, typename T1, typename T2>
void expect_ad(const F& f, const T1& x1, const T2& x2) {
  ad_tolerances tols;
  expect_ad(tols, f, x1, x2);
}

/**
 * Test that the specified ternary function produces autodiff values
 * and 1st-, 2nd-, and 3rd-order derivatives consistent with primitive
 * int and double inputs and finite differences at the specified
 * tolerances.
 *
 * @tparam F type of binary polymorphic functor to test
 * @tparam T1 type of double- or int-based first argument
 * @tparam T2 type of double- or int-based second argument
 * @tparam T3 type of double- or int-based third argument
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 first argument to test
 * @param x2 second argument to test
 * @param x3 third argument to test
 */
template <typename F, typename T1, typename T2, typename T3>
void expect_ad(const ad_tolerances& tols, const F& f, const T1& x1,
               const T2& x2, const T3& x3) {
  internal::expect_ad_vvv(tols, f, x1, x2, x3);
}

/**
 * Test that the specified ternary function produces autodiff values
 * and derivatives consistent with primitive int and double inputs and
 * finite differences at default tolerances.
 *
 * @tparam F type of binary polymorphic functor to test
 * @tparam T1 type of double- or int-based first argument
 * @tparam T2 type of double- or int-based second argument
 * @tparam T3 type of double- or int-based third argument
 * @param f functor to test
 * @param x1 first argument to test
 * @param x2 second argument to test
 * @param x3 third argument to test
 */
template <typename F, typename T1, typename T2, typename T3>
void expect_ad(const F& f, const T1& x1, const T2& x2, const T3& x3) {
  ad_tolerances tols;
  expect_ad(tols, f, x1, x2, x3);
}

/**
 * Promote to Complex is used by the `expect_ad_vectorized` framework
 *  to specify at compile time whether a function should check
 *  for complex support. By default the value is `Real`, meaning that
 *  the test does not promote the input to a complex double when it is running
 *  the test suite.
 */
enum class ScalarSupport { Real, RealAndComplex, ComplexOnly };

/**
 * Test that the specified vectorized polymorphic unary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for real numbers only.
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument (integer or double)
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 value to test
 */
template <
    ScalarSupport ComplexSupport = ScalarSupport::Real, typename F, typename T1,
    stan::require_t<
        stan::bool_constant<ComplexSupport == ScalarSupport::Real>>* = nullptr>
void expect_ad_vectorized(const ad_tolerances& tols, const F& f, const T1& x1) {
  using Scalar = std::conditional_t<std::is_integral<T1>::value, double, T1>;
  using matrix_t = Eigen::Matrix<Scalar, -1, -1>;
  using vector_t = Eigen::Matrix<Scalar, -1, 1>;
  using row_vector_t = Eigen::Matrix<Scalar, 1, -1>;
  using std::vector;
  typedef vector<Scalar> vector_dbl;
  typedef vector<vector<Scalar>> vector2_dbl;
  typedef vector<vector<vector<Scalar>>> vector3_dbl;

  expect_ad(tols, f, x1);
  expect_ad(tols, f, static_cast<Scalar>(x1));
  for (int i = 0; i < 2; ++i)
    expect_ad(tols, f, vector_t::Constant(i, x1).eval());
  for (int i = 0; i < 2; ++i)
    expect_ad(tols, f, row_vector_t::Constant(i, x1).eval());
  for (int i = 0; i < 2; ++i)
    expect_ad(tols, f, matrix_t::Constant(i, i, x1).eval());
  for (size_t i = 0; i < 2; ++i)
    expect_ad(tols, f, vector_dbl(i, x1));
  for (size_t i = 0; i < 2; ++i)
    expect_ad(tols, f, vector<vector_t>(i, vector_t::Constant(i, x1).eval()));
  for (size_t i = 0; i < 2; ++i)
    expect_ad(tols, f,
              vector<row_vector_t>(i, row_vector_t::Constant(i, x1).eval()));
  for (size_t i = 0; i < 2; ++i)
    expect_ad(tols, f,
              vector<matrix_t>(i, matrix_t::Constant(i, i, x1).eval()));
  for (int i = 0; i < 2; ++i)
    expect_ad(tols, f, vector2_dbl(i, vector_dbl(i, x1)));
  for (int i = 0; i < 2; ++i)
    expect_ad(tols, f, vector3_dbl(i, vector2_dbl(i, vector_dbl(i, x1))));
}

/**
 * Test that the specified vectorized polymorphic unary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when complex and real numbers
 * are supported.
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument (integer or double)
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 value to test
 */
template <ScalarSupport ComplexSupport, typename F, typename T1,
          stan::require_t<stan::bool_constant<
              ComplexSupport == ScalarSupport::RealAndComplex>>* = nullptr>
void expect_ad_vectorized(const ad_tolerances& tols, const F& f, const T1& x1) {
  using Scalar = std::conditional_t<std::is_integral<T1>::value, double, T1>;
  using matrix_t = Eigen::Matrix<Scalar, -1, -1>;
  using vector_t = Eigen::Matrix<Scalar, -1, 1>;
  using row_vector_t = Eigen::Matrix<Scalar, 1, -1>;
  using complex_t = std::complex<double>;
  using complex_matrix_t = Eigen::Matrix<complex_t, -1, -1>;
  using complex_vector_t = Eigen::Matrix<complex_t, -1, 1>;
  using complex_row_vector_t = Eigen::Matrix<complex_t, 1, -1>;
  using std::vector;
  using vector_dbl = vector<Scalar>;
  using vector2_dbl = vector<vector<Scalar>>;
  using vector3_dbl = vector<vector<vector<Scalar>>>;
  using vector_complex = vector<complex_t>;
  using vector2_complex = vector<vector_complex>;
  using vector3_complex = vector<vector2_complex>;

  expect_ad(tols, f, x1);
  expect_ad(tols, f, static_cast<Scalar>(x1));
  expect_ad(tols, f, std::complex<double>(static_cast<Scalar>(x1)));
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector_t::Constant(i, x1).eval());
    expect_ad(tols, f, complex_vector_t::Constant(i, x1).eval());
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, row_vector_t::Constant(i, x1).eval());
    expect_ad(tols, f, complex_row_vector_t::Constant(i, x1).eval());
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, matrix_t::Constant(i, i, x1).eval());
    expect_ad(tols, f, complex_matrix_t::Constant(i, i, x1).eval());
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector_dbl(i, x1));
    expect_ad(tols, f, vector_complex(i, x1));
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector<vector_t>(i, vector_t::Constant(i, x1).eval()));
    expect_ad(
        tols, f,
        vector<complex_vector_t>(i, complex_vector_t::Constant(i, x1).eval()));
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f,
              vector<row_vector_t>(i, row_vector_t::Constant(i, x1).eval()));
    expect_ad(tols, f,
              vector<complex_row_vector_t>(
                  i, complex_row_vector_t::Constant(i, x1).eval()));
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f,
              vector<matrix_t>(i, matrix_t::Constant(i, i, x1).eval()));
    expect_ad(tols, f,
              vector<complex_matrix_t>(
                  i, complex_matrix_t::Constant(i, i, x1).eval()));
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector2_dbl(i, vector_dbl(i, x1)));
    expect_ad(tols, f, vector2_complex(i, vector_complex(i, x1)));
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector3_dbl(i, vector2_dbl(i, vector_dbl(i, x1))));
    expect_ad(tols, f,
              vector3_complex(i, vector2_complex(i, vector_complex(i, x1))));
  }
}

/**
 * Test that the specified vectorized polymorphic unary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when only complex numbers are
 * supported.
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument (integer or double)
 * @param tols tolerances for test
 * @param f functor to test
 * @param x1 value to test
 */
template <ScalarSupport ComplexSupport, typename F, typename T1,
          stan::require_t<stan::bool_constant<
              ComplexSupport == ScalarSupport::ComplexOnly>>* = nullptr>
void expect_ad_vectorized(const ad_tolerances& tols, const F& f, const T1& x1) {
  using Scalar = std::conditional_t<std::is_integral<T1>::value, double, T1>;
  using complex_t = std::complex<double>;
  using complex_matrix_t = Eigen::Matrix<complex_t, -1, -1>;
  using complex_vector_t = Eigen::Matrix<complex_t, -1, 1>;
  using complex_row_vector_t = Eigen::Matrix<complex_t, 1, -1>;
  using std::vector;
  using vector_complex = vector<complex_t>;
  using vector2_complex = vector<vector_complex>;
  using vector3_complex = vector<vector2_complex>;

  expect_ad(tols, f, std::complex<double>(static_cast<Scalar>(x1)));
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, complex_vector_t::Constant(i, x1).eval());
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, complex_row_vector_t::Constant(i, x1).eval());
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, complex_matrix_t::Constant(i, i, x1).eval());
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector_complex(i, x1));
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(
        tols, f,
        vector<complex_vector_t>(i, complex_vector_t::Constant(i, x1).eval()));
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f,
              vector<complex_row_vector_t>(
                  i, complex_row_vector_t::Constant(i, x1).eval()));
  }
  for (size_t i = 0; i < 2; ++i) {
    expect_ad(tols, f,
              vector<complex_matrix_t>(
                  i, complex_matrix_t::Constant(i, i, x1).eval()));
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f, vector2_complex(i, vector_complex(i, x1)));
  }
  for (int i = 0; i < 2; ++i) {
    expect_ad(tols, f,
              vector3_complex(i, vector2_complex(i, vector_complex(i, x1))));
  }
}

/**
 * Test that the specified function has value and 1st-, 2nd-, and
 * 3rd-order derivatives consistent with primitive values and finite
 * differences using default tolerances.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both.
 * @tparam F type of function
 * @tparam T type of argument
 * @param f function to test
 * @param x argument to test
 */
template <ScalarSupport ComplexSupport = ScalarSupport::Real, typename F,
          typename T>
void expect_ad_vectorized(const F& f, const T& x) {
  ad_tolerances tols;
  expect_ad_vectorized<ComplexSupport>(tols, f, x);
}

/**
 * Implementation function for testing that binary functions with vector inputs
 * (both Eigen and std::vector types) return 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences of double inputs.
 *
 * @tparam F type of function
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 */
template <typename F, typename T1, typename T2,
          require_all_not_st_integral<T1, T2>* = nullptr>
void expect_ad_vectorized_binary_impl(const ad_tolerances& tols, const F& f,
                                      const T1& x, const T2& y) {
  std::vector<T1> nest_x{x, x};
  std::vector<T2> nest_y{y, y};
  std::vector<std::vector<T1>> nest_nest_x{nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y{nest_y, nest_y};
  expect_ad(tols, f, x, y);            // mat, mat
  expect_ad(tols, f, x, y[0]);         // mat, scal
  expect_ad(tols, f, x[0], y);         // scal, mat
  expect_ad(tols, f, nest_x, nest_y);  // nest<mat>, nest<mat>
  expect_ad(tols, f, nest_x, y[0]);    // nest<mat>, scal
  expect_ad(tols, f, x[0], nest_y);    // scal, nest<mat>
  expect_ad(tols, f, nest_nest_x,
            nest_nest_y);                 // nest<nest<mat>>, nest<nest<mat>>
  expect_ad(tols, f, nest_nest_x, y[0]);  // nest<nest<mat>, scal
  expect_ad(tols, f, x[0], nest_nest_y);  // scal, nest<nest<mat>>
}

/**
 * Implementation function for testing that ternary functions with vector inputs
 * (both Eigen and std::vector types) return 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences of double inputs.
 *
 * @tparam F type of function
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @tparam T3 type of third argument
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 * @param z argument to test
 */
template <typename F, typename T1, typename T2, typename T3>
void expect_ad_vectorized_ternary_impl(const ad_tolerances& tols, const F& f,
                                       const T1& x, const T2& y, const T3& z) {
  std::vector<T1> nest_x{x};
  std::vector<T2> nest_y{y};
  std::vector<T3> nest_z{z};
  std::vector<std::vector<T1>> nest_nest_x{nest_x};
  std::vector<std::vector<T2>> nest_nest_y{nest_y};
  std::vector<std::vector<T3>> nest_nest_z{nest_z};
  expect_ad(tols, f, x, y, z);
  expect_ad(tols, f, nest_nest_x, nest_nest_y, nest_nest_z);
  expect_ad(tols, f, nest_nest_x, nest_nest_y, z[0]);
  expect_ad(tols, f, nest_nest_x, y[0], z[0]);
  expect_ad(tols, f, nest_nest_x, y[0], nest_nest_z);
  expect_ad(tols, f, x[0], y[0], nest_nest_z);
  expect_ad(tols, f, x[0], nest_nest_y, nest_nest_z);
  expect_ad(tols, f, x[0], nest_nest_y, z[0]);
}

/**
 * Implementation function for testing that binary functions with vector inputs
 * (both Eigen and std::vector types) return 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences of double inputs.
 *
 * This is a specialisation for use when the first input is an integer type
 *
 * @tparam F type of function
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 */
template <typename F, typename T1, typename T2,
          require_st_integral<T1>* = nullptr>
void expect_ad_vectorized_binary_impl(const ad_tolerances& tols, const F& f,
                                      const T1& x, const T2& y) {
  auto f_bind
      = [&](const auto& x) { return [=](const auto& y) { return f(x, y); }; };
  std::vector<T1> nest_x{x, x};
  std::vector<T2> nest_y{y, y};
  std::vector<std::vector<T1>> nest_nest_x{nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y{nest_y, nest_y};
  expect_ad(tols, f_bind(x), y);
  expect_ad(tols, f_bind(nest_x), nest_y);
  expect_ad(tols, f_bind(nest_nest_x), nest_nest_y);
}

/**
 * Implementation function for testing that binary functions with vector inputs
 * (both Eigen and std::vector types) return 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences of double inputs.
 *
 * This is a specialisation for use when the second input is an integer type
 *
 * @tparam F type of function
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 */
template <typename F, typename T1, typename T2,
          require_st_integral<T2>* = nullptr>
void expect_ad_vectorized_binary_impl(const ad_tolerances& tols, const F& f,
                                      const T1& x, const T2& y) {
  auto f_bind
      = [&](const auto& y) { return [=](const auto& x) { return f(x, y); }; };
  std::vector<T1> nest_x{x, x};
  std::vector<T2> nest_y{y, y};
  std::vector<std::vector<T1>> nest_nest_x{nest_x, nest_x};
  std::vector<std::vector<T2>> nest_nest_y{nest_y, nest_y};
  expect_ad(tols, f_bind(y), x);
  expect_ad(tols, f_bind(nest_y), nest_x);
  expect_ad(tols, f_bind(nest_nest_y), nest_nest_x);
}

/**
 * Test that the specified vectorized polymorphic binary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument
 * @tparam T1 type of second argument
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param y value to test
 */
template <typename F, typename T1, typename T2,
          require_all_eigen_col_vector_t<T1, T2>* = nullptr>
void expect_ad_vectorized_binary(const ad_tolerances& tols, const F& f,
                                 const T1& x, const T2& y) {
  expect_ad_vectorized_binary_impl(tols, f, x, y);
  expect_ad_vectorized_binary_impl(tols, f, math::to_array_1d(x),
                                   math::to_array_1d(y));
}

/**
 * Test that the specified vectorized polymorphic ternary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @tparam T3 type of third argument
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param y value to test
 * @param z value to test
 */
template <typename F, typename T1, typename T2, typename T3,
          require_all_eigen_col_vector_t<T1, T2, T3>* = nullptr>
void expect_ad_vectorized_ternary(const ad_tolerances& tols, const F& f,
                                  const T1& x, const T2& y, const T3& z) {
  expect_ad_vectorized_ternary_impl(tols, f, x, y, z);
  expect_ad_vectorized_ternary_impl(tols, f, math::to_array_1d(x),
                                    math::to_array_1d(y), math::to_array_1d(z));
}

/**
 * Test that the specified vectorized polymorphic binary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument
 * @tparam T1 type of second argument
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param y value to test
 */
template <typename F, typename T1, typename T2,
          require_any_std_vector_t<T1, T2>* = nullptr>
void expect_ad_vectorized_binary(const ad_tolerances& tols, const F& f,
                                 const T1& x, const T2& y) {
  expect_ad_vectorized_binary_impl(tols, f, x, y);
}

/**
 * Test that the specified vectorized polymorphic ternary function
 * produces autodiff results consistent with values determined by
 * double and integer inputs and 1st-, 2nd-, and 3rd-order derivatives
 * consistent with finite differences of double inputs.
 *
 * @tparam F type of polymorphic, vectorized functor to test
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @tparam T3 type of third argument
 * @param tols tolerances for test
 * @param f functor to test
 * @param x value to test
 * @param y value to test
 * @param z value to test
 */
template <typename F, typename T1, typename T2, typename T3,
          require_any_std_vector_t<T1, T2, T3>* = nullptr>
void expect_ad_vectorized_ternary(const ad_tolerances& tols, const F& f,
                                  const T1& x, const T2& y, const T3& z) {
  expect_ad_vectorized_ternary_impl(tols, f, x, y, z);
}

/**
 * Test that the specified binary function has value and 1st-, 2nd-, and
 * 3rd-order derivatives consistent with primitive values and finite
 * differences using default tolerances.
 *
 * @tparam F type of function
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 */
template <typename F, typename T1, typename T2>
void expect_ad_vectorized_binary(const F& f, const T1& x, const T2& y) {
  ad_tolerances tols;
  expect_ad_vectorized_binary(tols, f, x, y);
}

/**
 * Test that the specified ternary function has value and 1st-, 2nd-, and
 * 3rd-order derivatives consistent with primitive values and finite
 * differences using default tolerances.
 *
 * @tparam F type of function
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @tparam T3 type of third argument
 * @param f function to test
 * @param x argument to test
 * @param y argument to test
 * @param z argument to test
 */
template <typename F, typename T1, typename T2, typename T3>
void expect_ad_vectorized_ternary(const F& f, const T1& x, const T2& y,
                                  const T3& z) {
  ad_tolerances tols;
  expect_ad_vectorized_ternary(tols, f, x, y, z);
}

/**
 * Test that the specified polymorphic unary function produces the
 * same results, exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function, when applied to the common
 * arguments.
 *
 * @tparam F type of polymorphic unary functor
 * @param f unary functor to test
 */
template <typename F>
void expect_common_unary(const F& f) {
  auto args = internal::common_args();
  for (double x1 : args)
    expect_ad(f, x1);
  auto int_args = internal::common_int_args();
  for (int x1 : int_args)
    expect_ad(f, x1);
}

/**
 * Test that the specified polymorphic unary function produces the
 * same results, exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function, when applied to all
 * common integer and double arguments excluding zero.
 *
 * @tparam F type of polymorphic binary functor
 * @param f functor to test
 */
template <typename F>
void expect_common_nonzero_unary(const F& f) {
  auto args = internal::common_nonzero_args();
  for (double x1 : args)
    expect_ad(f, x1);

  auto int_args = internal::common_nonzero_int_args();
  for (int x : int_args)
    expect_ad(f, x);
}

/**
 * Test that the specified polymorphic binary function produces the
 * same results, exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function, when applied to all pairs of
 * common integer and double argument combinations excluding zero.
 *
 * If the `disable_lhs_int` flag is set to `true` (it defaults to
 * `false`), then integers will not be considered as first arguments.
 * This is useful for testing assignment operators like `+=` and
 * division operators like `/` where integer and real arguments
 * produce different values.
 *
 * @tparam F type of polymorphic binary functor
 * @param f functor to test
 * @param disable_lhs_int if integer values should only be tested
 * for second argments
 */
template <typename F>
void expect_common_nonzero_binary(const F& f, bool disable_lhs_int = false) {
  auto args = internal::common_nonzero_args();
  auto int_args = internal::common_nonzero_int_args();
  for (double x1 : args)
    for (double x2 : args) {
      expect_ad(f, x1, x2);
    }
  for (double x1 : args)
    for (int x2 : int_args) {
      expect_ad(f, x1, x2);
    }

  if (disable_lhs_int)
    return;

  for (int x1 : int_args)
    for (double x2 : args) {
      expect_ad(f, x1, x2);
    }
  for (int x1 : int_args)
    for (int x2 : int_args) {
      expect_ad(f, x1, x2);
    }
}

/**
 * Test that the specified polymorphic binary function produces the
 * same results, exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function, when applied to all pairs of
 * common integer and double argument combinations.
 *
 * If the `disable_lhs_int` flag is set to `true` (it defaults to
 * `false`), then integers will not be considered as first arguments.
 * This is useful for testing assignment operators like `+=` and
 * division operators like `/` where integer and real arguments
 * produce different values.
 *
 * @tparam F type of polymorphic binary functor
 * @param f functor to test
 * @param disable_lhs_int if integer values should only be tested
 * for second argments
 */
template <typename F>
void expect_common_binary(const F& f, bool disable_lhs_int = false) {
  auto args = internal::common_args();
  auto int_args = internal::common_int_args();
  for (double x1 : args)
    for (double x2 : args) {
      expect_ad(f, x1, x2);
    }
  for (double x1 : args)
    for (int x2 : int_args) {
      expect_ad(f, x1, x2);
    }
  if (disable_lhs_int)
    return;
  for (int x1 : int_args)
    for (double x2 : args) {
      expect_ad(f, x1, x2);
    }
  for (int x1 : int_args)
    for (int x2 : int_args) {
      expect_ad(f, x1, x2);
    }
}

std::vector<double> common_complex_parts_re() {
  return {-4, -2.5, -1.5, -0.3, -0.1, 0.1, 1.3, 2.1, 3.9};
}

std::vector<double> common_complex_parts_im() {
  return {-4, -2.5, -1.5, -0.3, -0.0, 0.0, 1.3, 2.1, 3.9};
}

std::vector<std::complex<double>> common_complex() {
  std::vector<std::complex<double>> zs;
  auto complex_re = common_complex_parts_re();
  auto complex_im = common_complex_parts_im();
  for (int i = 0; i < complex_re.size(); ++i) {
    zs.emplace_back(complex_re[i], complex_im[i]);
  }
  return zs;
}

template <typename F>
void expect_complex_common(const F& f) {
  auto zs = common_complex();
  for (auto z : zs) {
    expect_ad(f, z);
  }
}

template <typename F>
void expect_complex_common_binary(const F& f) {
  auto xs = common_complex_parts_re();
  auto zs = common_complex();
  // complex, complex
  for (auto z1 : zs) {
    for (auto z2 : zs) {
      expect_ad(f, z1, z2);
    }
  }
  // complex, real
  for (auto z1 : zs) {
    for (auto x2 : xs) {
      expect_ad(f, z1, x2);
    }
  }
  // real, complex
  for (auto x1 : xs) {
    for (auto z2 : zs) {
      expect_ad(f, x1, z2);
    }
  }
}

template <typename T, typename F>
void expect_complex_compare(const F& f, const std::complex<double>& z1,
                            const std::complex<double>& z2) {
  using c_t = std::complex<T>;
  c_t cz1{z1};
  c_t cz2{z2};
  T z1r{z1.real()};
  T z2r{z2.real()};

  EXPECT_EQ(f(z1, z2), f(cz1, cz2));
  EXPECT_EQ(f(z1, z2), f(cz1, z2));
  EXPECT_EQ(f(z1, z2), f(z1, cz2));

  EXPECT_EQ(f(z1.real(), z2), f(z1r, cz2));
  EXPECT_EQ(f(z1.real(), z2), f(z1r, z2));

  EXPECT_EQ(f(z1, z2.real()), f(cz1, z2r));
  EXPECT_EQ(f(z1, z2.real()), f(z1, z2r));
}

template <typename F>
void expect_complex_comparison(const F& f, const std::complex<double>& z1,
                               const std::complex<double>& z2) {
  using stan::math::fvar;
  using stan::math::var;
  using std::complex;
  expect_complex_compare<double>(f, z1, z2);
  expect_complex_compare<var>(f, z1, z2);
  expect_complex_compare<fvar<double>>(f, z1, z2);
  expect_complex_compare<fvar<fvar<double>>>(f, z1, z2);
  expect_complex_compare<fvar<var>>(f, z1, z2);
  expect_complex_compare<fvar<fvar<var>>>(f, z1, z2);
}

/**
 * Test the specified comparison operation provides results matching
 * those for the double version for all the common complex numbers.
 *
 * @tparam F type of function to test
 * @param f function to test
 */
template <typename F>
void expect_complex_common_comparison(const F& f) {
  for (auto z1 : common_complex()) {
    for (auto z2 : common_complex()) {
      expect_complex_comparison(f, z1, z2);
    }
  }
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function when applied to all common
 * arguments.  Uses default tolerances.
 *
 * <p>The function must be defined from scalars to scalars and from
 * containers to containers, always producing the same output type as
 * input type.  The value for containers must be the same as applying
 * the scalar function elementwise.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for whenever complex numbers are not
 * supported.
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <
    ScalarSupport ComplexSupport = ScalarSupport::Real, typename F,
    require_t<bool_constant<ComplexSupport == ScalarSupport::Real>>* = nullptr>
void expect_common_unary_vectorized(const F& f) {
  ad_tolerances tols;
  auto args = internal::common_args();
  for (double x1 : args)
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
  auto int_args = internal::common_int_args();
  for (int x1 : args)
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function when applied to all common
 * arguments.  Uses default tolerances.
 *
 * <p>The function must be defined from scalars to scalars and from
 * containers to containers, always producing the same output type as
 * input type.  The value for containers must be the same as applying
 * the scalar function elementwise.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when complex and real numbers
 * are supported.
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <ScalarSupport ComplexSupport, typename F,
          require_t<bool_constant<ComplexSupport
                                  == ScalarSupport::RealAndComplex>>* = nullptr>
void expect_common_unary_vectorized(const F& f) {
  ad_tolerances tols;
  auto args = internal::common_args();
  for (double x1 : args)
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
  auto int_args = internal::common_int_args();
  for (int x1 : args)
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
  for (auto x1 : common_complex())
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has 1st-, 2nd-, and 3rd-order
 * derivatives consistent with finite differences as returned by the
 * primitive version of the function when applied to all common
 * arguments.  Uses default tolerances.
 *
 * <p>The function must be defined from scalars to scalars and from
 * containers to containers, always producing the same output type as
 * input type.  The value for containers must be the same as applying
 * the scalar function elementwise.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when only complex numbers are
 * supported.
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <ScalarSupport ComplexSupport, typename F,
          require_t<bool_constant<ComplexSupport
                                  == ScalarSupport::ComplexOnly>>* = nullptr>
void expect_common_unary_vectorized(const F& f) {
  ad_tolerances tols;
  for (auto x1 : common_complex())
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
}

template <ScalarSupport ComplexSupport = ScalarSupport::Real, typename F>
void expect_unary_vectorized(const ad_tolerances& tols, const F& f) {}

/**
 * Test that the specified vectorized unary function has value and
 * derivative behavior matching the primitive instantiation with finite
 * differences.  Tests both scalar and container behavior.  Integer
 * arguments will be preserved through to function calls.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both.
 * @tparam F type of function
 * @tparam Ts types of arguments
 * @param tols test relative tolerances
 * @param f function to test
 * @param xs arguments to test
 */
template <ScalarSupport ComplexSupport = ScalarSupport::Real, typename F,
          typename T, typename... Ts>
void expect_unary_vectorized(const ad_tolerances& tols, const F& f, T x,
                             Ts... xs) {
  expect_ad_vectorized<ComplexSupport>(tols, f, x);
  expect_unary_vectorized<ComplexSupport>(tols, f, xs...);
}

/**
 * Test that the specified unary function produces derivatives and
 * values for the specified values that are consistent with primitive
 * values and finite differences.  Tests both scalars and containers.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both.
 * @tparam F type of function to test
 * @tparam Ts type of remaining arguments to test
 * @param f function to test
 * @param xs arguments to test
 */
template <ScalarSupport ComplexSupport = ScalarSupport::Real, typename F,
          require_not_same_t<F, ad_tolerances>* = nullptr, typename... Ts>
void expect_unary_vectorized(const F& f, Ts... xs) {
  ad_tolerances tols;  // default tolerances
  expect_unary_vectorized<ComplexSupport>(tols, f, xs...);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has derivatives consistent with finite
 * differences as returned by the primitive version of the function
 * when applied to all common non-zero integer and double arguments.
 * This includes tests for standard vector and Eigen vector containers.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when complex numbers are not
 * supported.
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <ScalarSupport ComplexSupport = ScalarSupport::Real, typename F,
          stan::require_t<stan::bool_constant<
              ComplexSupport == ScalarSupport::Real>>* = nullptr>
void expect_common_nonzero_unary_vectorized(const F& f) {
  ad_tolerances tols;
  for (double x : internal::common_nonzero_args())
    stan::test::expect_unary_vectorized<ComplexSupport>(tols, f, x);
  for (auto x : internal::common_nonzero_int_args())
    stan::test::expect_unary_vectorized<ComplexSupport>(tols, f, x);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has derivatives consistent with finite
 * differences as returned by the primitive version of the function
 * when applied to all common non-zero integer and double arguments.
 * This includes tests for standard vector and Eigen vector containers.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when complex and real numbers
 * are supported.
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <ScalarSupport ComplexSupport, typename F,
          stan::require_t<stan::bool_constant<
              ComplexSupport == ScalarSupport::RealAndComplex>>* = nullptr>
void expect_common_nonzero_unary_vectorized(const F& f) {
  ad_tolerances tols;
  for (double x : internal::common_nonzero_args())
    stan::test::expect_unary_vectorized<ComplexSupport>(tols, f, x);
  for (int x : internal::common_nonzero_int_args())
    stan::test::expect_unary_vectorized<ComplexSupport>(tols, f, x);
  for (auto x1 : common_complex())
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has derivatives consistent with finite
 * differences as returned by the primitive version of the function
 * when applied to all common non-zero integer and double arguments.
 * This includes tests for standard vector and Eigen vector containers.
 *
 * @tparam ScalarSupport whether the input supports real numbers, complex
 * numbers, or both. This specialization is for when only complex numbers are
 * supported.
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <ScalarSupport ComplexSupport, typename F,
          stan::require_t<stan::bool_constant<
              ComplexSupport == ScalarSupport::ComplexOnly>>* = nullptr>
void expect_common_nonzero_unary_vectorized(const F& f) {
  ad_tolerances tols;
  for (auto x1 : common_complex())
    stan::test::expect_ad_vectorized<ComplexSupport>(tols, f, x1);
}

/**
 * For all pairs of common arguments, test that primitive and all
 * autodiff types return the same value. Tests integer and double
 * arguments.
 *
 * @tparam F type of polymorphic binary functor
 * @param f functor to test
 */
template <typename F>
void expect_common_comparison(const F& f) {
  auto args = internal::common_args();
  auto int_args = internal::common_int_args();
  for (double x1 : args)
    for (double x2 : args)
      internal::expect_comparison(f, x1, x2);
  for (int x1 : int_args)
    for (double x2 : args)
      internal::expect_comparison(f, x1, x2);
  for (double x1 : args)
    for (int x2 : int_args)
      internal::expect_comparison(f, x1, x2);
  for (int x1 : int_args)
    for (int x2 : int_args)
      internal::expect_comparison(f, x1, x2);
}

/**
 * Test that the two specified functions either both throw or have the
 * same return value for the specified argument.  If the argument is
 * an integer, it will be passed through to the functions as such.
 *
 * @tparam F1 type of first function
 * @tparam F2 type of second function
 * @tparam T type of argument
 * @param f1 first function to test
 * @param f2 second function to test
 * @param x argument to test
 */
template <typename F1, typename F2, typename T>
void expect_match_prim(const F1& f1, const F2& f2, const T& x) {
  try {
    auto y1 = f1(x);
    try {
      auto y2 = f2(x);
      // neither throw, so expect values to be the same
      expect_near_rel("expect_match_prim", y1, y2);
      SUCCEED() << "expect_match_prim: f1 and f2 return same values";
    } catch (...) {
      FAIL() << "expect_match_prim: f2 throws but f1 does not";
    }
  } catch (...) {
    try {
      f2(x);
      FAIL() << "expect_match_prim: f1 throws but f2 does not";
    } catch (...) {
      SUCCEED() << "expect_match_prim: f1 and f2 both throw";
      return;
    }
  }
}

/**
 * Test that the specified pair of functions either both throw or
 * return the same value for common integer and double inputs
 * including NaN and infinities.  Tests integer and double arguments.
 *
 * @tparam F1 type of first function
 * @tparam F2 type of second function
 * @param f1 first function to test
 * @param f2 second function to test
 */
template <typename F1, typename F2>
void expect_common_prim(const F1& f1, const F2& f2) {
  for (double x : internal::common_args())
    expect_match_prim(f1, f2, x);
  for (int x : internal::common_int_args())
    expect_match_prim(f1, f2, x);
}

/**
 * Return sequence of covariance matrices ranging in dimension from
 * the specified minimum to maximum, with specified autocorrelation.
 * The entries `Sigma[i, j]` in the result are `pow(rho, fabs(i -
 * j))`.
 *
 * @param N_min minimum covariance matrix dimension
 * @param N_max maximum covariance matrix dimension
 * @param rho correlation (between -1 and 1)
 * @return sequence of covariance matrices between specified sizes
 * with specified autocorrelation
 */
std::vector<Eigen::MatrixXd> ar_test_cov_matrices(int N_min, int N_max,
                                                  double rho) {
  std::vector<Eigen::MatrixXd> ys;
  for (int n = N_min; n <= N_max; ++n) {
    Eigen::MatrixXd y(n, n);
    for (int i = 0; i < n; ++i) {
      y(i, i) = 1;
      for (int j = 0; j < i; ++j) {
        y(i, j) = std::pow(rho, std::fabs(i - j));
        y(j, i) = y(i, j);
      }
    }
    ys.push_back(y);
  }
  return ys;
}

/**
 * Return Eigen vector with elements given by the specified standard
 * vector.
 *
 * @param x standard vector input
 * @return copy as Eigen vector
 */
Eigen::VectorXd to_vector(const std::vector<double>& x) {
  Eigen::VectorXd y(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    y(i) = x[i];
  return y;
}

/**
 * Copy the specified Eigen matrix, vector, or row vector to a
 * vector.
 *
 * @tparam R row specification for matrix
 * @tparam C column specification for matrix
 * @param x matrix
 * @return copy as vector
 */
template <int R, int C>
Eigen::VectorXd to_vector(const Eigen::Matrix<double, R, C>& x) {
  Eigen::VectorXd y(x.size());
  for (int i = 0; i < x.size(); ++i)
    y(i) = x(i);
  return y;
}

/**
 * Return Eigen row vector with elements given by the specified
 * standard vector.
 *
 * @param x standard vector input
 * @return copy as Eigen row vector
 */
Eigen::RowVectorXd to_row_vector(const std::vector<double>& x) {
  Eigen::RowVectorXd y(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    y(i) = x[i];
  return y;
}

/**
 * Copy the specified Eigen matrix, vector, or row vector to a
 * row vector.
 *
 * @tparam R row specification for matrix
 * @tparam C column specification for matrix
 * @param x matrix
 * @return copy as row vector
 */
template <int R, int C>
Eigen::VectorXd to_row_vector(const Eigen::Matrix<double, R, C>& x) {
  Eigen::RowVectorXd y(x.size());
  for (int i = 0; i < x.size(); ++i)
    y(i) = x(i);
  return y;
}

/**
 * Return square test matrices (not symmetric) of dimensionality
 * within the specified range (inclusive).
 *
 * @param min minimum matrix dimensionality to include
 * @param max maximum matrix dimensionality to include
 * @return square matrices within given dimensionality range (inclusive)
 */
std::vector<Eigen::MatrixXd> square_test_matrices(int low, int high) {
  std::vector<Eigen::MatrixXd> xs;
  Eigen::MatrixXd a00(0, 0);
  if (0 >= low && 0 <= high)
    xs.push_back(a00);

  Eigen::MatrixXd a11(1, 1);
  a11 << -1.3;
  if (1 >= low && 1 <= high)
    xs.push_back(a11);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 0.7;
  if (2 >= low && 2 <= high)
    xs.push_back(a22);

  Eigen::MatrixXd a33(3, 3);
  a33 << 3, -5, 7, -7.2, 9.1, -6.3, 7, 12, -3;
  if (3 >= low && 3 <= high)
    xs.push_back(a33);

  return xs;
}

}  // namespace test
}  // namespace stan
#endif
