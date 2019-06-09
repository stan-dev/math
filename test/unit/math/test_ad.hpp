#ifndef TEST_UNIT_MATH_TEST_AD_HPP
#define TEST_UNIT_MATH_TEST_AD_HPP

#include <stan/math/mix/mat.hpp>
#include <test/unit/math/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>

namespace stan {
namespace test {

/**
 * Return true if the specified value is finite.
 *
 * @param x value to test
 * @return true if value is finite
 */
bool is_finite(double x) {
  return !stan::math::is_inf(x) && !stan::math::is_nan(x);
}

/**
 * Return true if all of the elements in the container are finite
 *
 * @tparam T scalar type
 * @tparam R row type
 * @tparam C col type
 * @param x container to test
 * @return true if all container values are finite
 */
template <typename T, int R, int C>
bool is_finite(const Eigen::Matrix<T, R, C>& x) {
  for (int i = 0; i < x.size(); ++i)
    if (!is_finite(x(i)))
      return false;
  return true;
}

/**
 * Return true if all of the elements in the container are finite
 *
 * @tparam T contained type
 * @param x container to test
 * @return true if all container values are finite
 */
template <typename T>
bool is_finite(const std::vector<T>& x) {
  for (size_t i = 0; i < x.size(); ++i)
    if (!is_finite(x[i]))
      return false;
  return true;
}

/**
 * Test that the specified values are within the specified tolerance
 * on relative error, and if not, fail the embedded google test.
 *
 * <p>Relative error is defined to be the error `u - v` rescaled by the
 * average absolute value,
 * `rel_err(u, v) = (u - v) / (0.5 * (abs(u) * + abs(v))).`
 *
 * <p>If both `u` and `v` are zero, the error is defined to be zero.
 * If only one of `u` or `v` is zero, the above definition always
 * yields 1.  This is not a useful test, so we instead test absolute
 * error in the case that one of the values is zero.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param msg failure message for reporting
 * @param x1 first argument
 * @param x2 second argument
 * @param tol relative tolerance
 */
template <typename T1, typename T2>
void expect_near_relative(const std::string& msg, const T1& x1, const T2& x2,
                          double tol = 1e-8) {
  using stan::math::fabs;
  if (x1 == 0 && x2 == 0)
    return;
  if (x1 == 0 || x2 == 0) {
    EXPECT_NEAR(x1, x2, tol) << "expect_near_relative(" << x1 << ", " << x2
                             << ", tolerance = " << tol << ")"
                             << "    in: " << msg << std::endl;
    return;
  }
  auto avg = 0.5 * (fabs(x1) + fabs(x2));
  auto relative_diff = (x1 - x2) / avg;
  EXPECT_NEAR(0, relative_diff, tol)
      << "expect_near_relative(" << x1 << ", " << x2 << ", tolerance = " << tol
      << ")"
      << "; relative diff = " << relative_diff << std::endl
      << "    in: " << msg << std::endl;
}

/**
 * Test that scalars x1 and x2 are within the specified relative
 * tolerance, with identity behavior for infinite and NaN values.
 * Relative tolerance is defined in the documentation for
 * `expect_near_relative`.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param msg message indicating what is being tested to print in case
 * of error
 * @param x1 first argument to test
 * @param x2 second argument to test
 * @param tol relative tolerance
 */
template <typename T1, typename T2>
void expect_near(const std::string& msg, const T1& x1, const T2& x2,
                 double tol = 1e-8) {
  if (stan::math::is_nan(x1) || stan::math::is_nan(x2))
    EXPECT_TRUE(stan::math::is_nan(x1) && stan::math::is_nan(x2))
        << "expect_near(" << x1 << ", " << x2 << ")" << std::endl
        << msg << std::endl;
  else if (stan::math::is_inf(x1) || stan::math::is_inf(x2))
    EXPECT_EQ(x1, x2) << "expect_near(" << x1 << ", " << x2 << ")" << std::endl
                      << msg << std::endl;
  else
    expect_near_relative(msg, x1, x2, tol);
}

/**
 * Tests that matrices (or vectors) x1 and x2 are same size and
 * have near values up to specified relative tolerance, accounting for
 * infinite, not-a-number, and zero values as specified in
 * `expect_near_relative`.
 *
 * @tparam T type of scalars in matrix
 * @tparam R row specification of matrix
 * @tparam C column specification of the matrix
 * @param msg failure message
 * @param x1 first matrix to test
 * @param x2 second matrix to test
 * @param tol relative tolerance
 */
template <typename T, int R, int C>
void expect_near(const std::string& msg, const Eigen::Matrix<T, R, C>& x1,
                 const Eigen::Matrix<T, R, C>& x2, double tol = 1e-8) {
  EXPECT_EQ(x1.rows(), x2.rows()) << "expect_near rows expect_eq(" << x1.rows()
                                  << ", " << x2.rows() << ")" << std::endl
                                  << msg << std::endl;
  EXPECT_EQ(x1.cols(), x2.cols()) << "expect_near cols expect_eq(" << x1.rows()
                                  << ", " << x2.rows() << ")" << std::endl
                                  << msg << std::endl;
  std::string msg2 = "expect_near elt x1(i) = x2(i)\n" + msg;
  for (int i = 0; i < x1.size(); ++i)
    expect_near(msg2, x1(i), x2(i), tol);
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
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_gradient(const F& f, const Eigen::VectorXd& x, double fx,
                   bool test_derivs = true) {
  Eigen::VectorXd grad_ad;
  double fx_ad = fx;
  stan::math::gradient<F>(f, x, fx_ad, grad_ad);
  expect_near("test_gradient fx = fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  stan::math::finite_diff_gradient(f, x, fx_fd, grad_fd);
  expect_near("test gradient grad_fd == grad_ad", grad_fd, grad_ad, 1e-4);
}

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
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_gradient_fvar(const F& f, const Eigen::VectorXd& x, double fx,
                        bool test_derivs = true) {
  Eigen::VectorXd grad_ad;
  double fx_ad = fx;
  stan::math::gradient<double, F>(f, x, fx_ad, grad_ad);
  expect_near("gradient_fvar fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  stan::math::finite_diff_gradient(f, x, fx_fd, grad_fd);
  expect_near("gradient_fvar grad_fd == grad_ad", grad_fd, grad_ad, 1e-4);
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
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_hessian_fvar(const F& f, const Eigen::VectorXd& x, double fx,
                       bool test_derivs = true) {
  double fx_ad;
  Eigen::VectorXd grad_ad;
  Eigen::MatrixXd H_ad;
  stan::math::hessian<double, F>(f, x, fx_ad, grad_ad, H_ad);
  expect_near("hessian_fvar fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  stan::math::finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
  expect_near("hessian fvar grad_fd == grad_ad", grad_fd, grad_ad, 1e-4);
  expect_near("hessian fvar H_fd = H_ad", H_fd, H_ad, 1e-3);
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
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_hessian(const F& f, const Eigen::VectorXd& x, double fx,
                  bool test_derivs = true) {
  double fx_ad;
  Eigen::VectorXd grad_ad;
  Eigen::MatrixXd H_ad;
  stan::math::hessian<F>(f, x, fx_ad, grad_ad, H_ad);
  expect_near("hessian fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  stan::math::finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
  expect_near("hessian grad_fd = grad_ad", grad_fd, grad_ad, 1e-4);
  expect_near("hessian grad_fd H_fd == H_ad", H_fd, H_ad, 1e-3);
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
 * @param f functor to test
 * @param x value to test
 * @param fx expected value
 * @param test_derivs `true` if derivatives should be tested
 */
template <typename F>
void test_grad_hessian(const F& f, const Eigen::VectorXd& x, double fx,
                       bool test_derivs = true) {
  double fx_ad;
  Eigen::MatrixXd H_ad;
  std::vector<Eigen::MatrixXd> grad_H_ad;
  stan::math::grad_hessian(f, x, fx_ad, H_ad, grad_H_ad);
  expect_near("grad_hessian fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::MatrixXd H_fd;
  std::vector<Eigen::MatrixXd> grad_H_fd;
  stan::math::finite_diff_grad_hessian(f, x, fx_fd, H_fd, grad_H_fd);
  expect_near("grad hessian H_fd == H_ad", H_fd, H_ad, 1e-3);
  EXPECT_EQ(x.size(), grad_H_fd.size());
  for (size_t i = 0; i < grad_H_fd.size(); ++i)
    expect_near("grad hessian grad_H_fd[i] == grad_H_ad[i]", grad_H_fd[i],
                grad_H_ad[i], 1e-2);
}

/**
 * For the specified functor and argument, test that automatic
 * differentiation provides the value as the double-based version and
 * the same derivatives as finite differences over the double
 * version.  It also tests that an autodiff version throws an
 * exception if and only if the double-based version throws an
 * exception.  It also ensures consistency of infinite, not-a-number,
 * and zero results as defined by the function `expect_near`.
 *
 * <p>The functor to test must define `T operator()(const
 * Matrix<T, -1, 1>& x) const;` for relevant `double` and autodiff types.
 *
 * <p>All results are tested for relative tolerance at thresholds
 * `1e-8` for values, `1e-4` for gradients (1st order derivatives),
 * `1e-3` for Hessians (2nd order derivatives), and `1e-2` for
 * gradients of Hessians (3rd order derivatives).  The reason for such
 * seemingly lax tests is that finite differences can be highly
 * unstable.
 *
 * @tparam G type of polymorphic functor
 * @param g polymorphic functor from vectors to scalars
 * @param x argument to test
 */
template <typename G>
void expect_ad_derivatives(const G& g, const Eigen::VectorXd& x) {
  double gx = g(x);
  test_gradient(g, x, gx);
  test_gradient_fvar(g, x, gx);
  test_hessian(g, x, gx);
  test_hessian_fvar(g, x, gx);
  test_grad_hessian(g, x, gx);
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
 * @param x values to test
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

template <typename F>
void expect_all_throw(const F& f, const Eigen::VectorXd& x) {
  using stan::math::fvar;
  using stan::math::var;
  expect_throw<var>(f, x, "var");
  expect_throw<fvar<double>>(f, x, "fvar<double>");
  expect_throw<fvar<fvar<double>>>(f, x, "fvar<fvar<double>>");
  expect_throw<fvar<var>>(f, x, "fvar<var>");
  expect_throw<fvar<fvar<var>>>(f, x, "fvar<fvar<var>>");
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
 * @param f functor to evaluate
 * @param h serialized functor taking an index and returning the function
 * that returns that component of the original output
 * original output
 * @param x serialized input
 * @param xs sequence of arguments with double-based scalars
 */
template <typename F, typename H, typename... Ts>
void expect_ad_helper(const F& f, const H& h, const Eigen::VectorXd& x,
                      Ts... xs) {
  size_t result_size;
  try {
    auto y = f(xs...);
    result_size = serialize<double>(y).size();
  } catch (...) {
    expect_all_throw(h(0), x);
    return;
  }
  for (size_t i = 0; i < result_size; ++i)
    expect_ad_derivatives(h(i), x);
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
 * @param x1 first argument
 * @param x2 second argument
 */
template <typename F, typename T1, typename T2>
void expect_ad_vv(const F& f, const T1& x1, const T2& x2) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x1ds = ds.read(x1);
      auto x2ds = ds.read(x2);
      return serialize_return(f(x1ds, x2ds))[i];
    };
  };
  expect_ad_helper(f, h, serialize_args(x1, x2), x1, x2);
}

/**
 * Test that the specified binary functor and arguments produce for
 * every autodiff type the same value as the double-based version and
 * the same derivatives as finite differences when the first argument
 * is an autodiff variable and the second double-based.
 *
 * @tparam F type of functor to test
 * @tparam T1 type of first argument with double-based scalar
 * @tparam T2 type of second argument with double-based scalar
 * @param f functor to test
 * @param x1 first argument
 * @param x2 second argument
 */
template <typename F, typename T1, typename T2>
void expect_ad_vd(const F& f, const T1& x1, const T2& x2) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x1ds = ds.read(x1);
      return serialize_return(f(x1ds, x2))[i];
    };
  };
  Eigen::VectorXd x = serialize_args(x1);
  expect_ad_helper(f, h, serialize_args(x1), x1, x2);
}

/**
 * Test that the specified binary functor and arguments produce for
 * every autodiff type the same value as the double-based version and
 * the same derivatives as finite differences when the second argument
 * is an autodiff variable and the first is double-based.
 *
 * @tparam F type of functor to test
 * @tparam T1 type of first argument with double-based scalar
 * @tparam T2 type of second argument with double-based scalar
 * @param f functor to test
 * @param x1 first argument
 * @param x2 second argument
 */
template <typename F, typename T1, typename T2>
void expect_ad_dv(const F& f, const T1& x1, const T2& x2) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x2ds = ds.read(x2);
      return serialize_return(f(x1, x2ds))[i];
    };
  };
  expect_ad_helper(f, h, serialize_args(x2), x1, x2);
}

/**
 * Test that the specified unary functor and arguments produce for
 * every autodiff type the same value as the double-based version and
 * the same derivatives as finite differences when the first (and
 * only) argument is an autodiff variable.
 *
 * @tparam F type of functor to test
 * @tparam T1 type of first argument with double-based scalar
 * @param f functor to test
 * @param x1 first argument
 */
template <typename F, typename T1>
void expect_ad_v(const F& f, const T1& x1) {
  auto h = [&](const int& i) {
    return [&](const auto& v) {
      auto ds = to_deserializer(v);
      auto x1ds = ds.read(x1);
      return serialize_return(f(x1ds))[i];
    };
  };
  expect_ad_helper(f, h, serialize_args(x1), x1);
}

/**
 * Test that the specified polymorphic unary functor produces autodiff
 * results consistent with values determined by double inputs and
 * derivatives consistent with finite differences of double inputs.
 *
 * <p>Tests condition where argument is an autodiff variable.  Tests
 * autodiff levels `rev`, `fvar<double>`, `fvar<fvar<double>>`,
 * `fvar<rev>`, and `fvar<fvar<rev>>`.
 *
 * <p>Invokes Google test framework to raise error if test fails.
 *
 * @tparam F type of functor to test
 * @tparam T type of argument
 */
template <typename F, typename T>
void expect_ad(const F& f, const T& x) {
  expect_ad_v(f, x);
}

/**
 * Test that the specified polymorphic binary functor produces autodiff
 * results consistent with values determined by double inputs and
 * derivatives consistent with finite differences of double inputs.
 *
 * <p>Comparison operations (`operator==`, `operator!=`, etc.) are
 * step functions when their inputs are equivalent, so their
 * derivatives are undefined and should not be tested via finite
 * differences.  The tests for derivatives can be turned off by
 * setting the final argument `is_comparison` to `true`; it takes a
 * default value of `false`.
 *
 * <p>Tests all three possible instantiations of autodiff variables:
 * first argument only, second argument only, and both arguments.
 * Tests autodiff levels `rev`, `fvar<double>`, `fvar<fvar<double>>`,
 * `fvar<rev>`, and `fvar<fvar<rev>>`.
 *
 * <p>Invokes Google test framework to raise error if test fails.
 *
 * @tparam F type of binary polymorphic functor to test
 * @tparam T1 type of double- or int-based first argument
 * @tparam T2 type of double- or int-based second argument
 * @param f functor to test
 * param x1 first argument to test
 * @param x2 second argument to test
 */
template <typename F, typename T1, typename T2>
void expect_ad(const F& f, const T1& x1, const T2& x2) {
  expect_ad_vv(f, x1, x2);
  expect_ad_vd(f, x1, x2);
  expect_ad_dv(f, x1, x2);
}

/**
 * Test that the specified vectorized polymoprhic unary function
 * produces autodiff results consistent with values determined by
 * double in puts and derivatives consistent with finite differences
 * of double inputs.
 *
 * <p>Tests all three possible instantiations of autodiff variables:
 * first argument only, second argument only, and both arguments.
 * Tests autodiff levels `rev`, `fvar<double>`, `fvar<fvar<double>>`,
 * `fvar<rev>`, and `fvar<fvar<rev>>`.
 *
 * <p>Tests all vectorizations of the second argument, including
 * primitive (`double` or `int`), `std::vector<double>` and all of the
 * Eigen options, `Eigen::VectorXd`, `Eigen::RowVectorXd`, and
 * `Eigen::MatrixXd`.   The vectorization tests are carried out by
 * repeating the input multiple times.

 * <p>Invokes Google test framework to raise error if test fails.
 *
 * @tparam F type of poymorphic, vectorized functor to test
 * @tparam T1 type of first argument (integer or double)
 * @param f functor to test
 * @param x1 value to test
 */
template <typename F, typename T1>
void expect_ad_vectorized(const F& f, const T1& x1) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using std::vector;
  typedef vector<double> vector_dbl;
  typedef vector<vector<double>> vector2_dbl;
  typedef vector<vector<vector<double>>> vector3_dbl;

  expect_ad(f, x1);
  expect_ad(f, static_cast<double>(x1));
  for (int i = 0; i < 4; ++i)
    expect_ad(f, VectorXd::Constant(i, x1).eval());
  for (int i = 0; i < 4; ++i)
    expect_ad(f, RowVectorXd::Constant(i, x1).eval());
  for (int i = 0; i < 4; ++i)
    expect_ad(f, MatrixXd::Constant(i, i, x1).eval());
  for (size_t i = 0; i < 4; ++i)
    expect_ad(f, vector_dbl(i, x1));
  for (size_t i = 0; i < 4; ++i)
    expect_ad(f, vector<VectorXd>(i, VectorXd::Constant(i, x1).eval()));
  for (size_t i = 0; i < 4; ++i)
    expect_ad(f, vector<RowVectorXd>(i, RowVectorXd::Constant(i, x1).eval()));
  for (size_t i = 0; i < 3; ++i)
    expect_ad(f, vector<MatrixXd>(i, MatrixXd::Constant(i, i, x1).eval()));
  for (int i = 0; i < 3; ++i)
    expect_ad(f, vector2_dbl(i, vector_dbl(i, x1)));
  for (int i = 0; i < 3; ++i)
    expect_ad(f, vector3_dbl(i, vector2_dbl(i, vector_dbl(i, x1))));
}

/**
 * Return a sequence of common non-zero arguments.  This includes
 * positive, negative, positive infinite, negative infinity, and
 * not-a-number values, but does not include zero.
 *
 * @return non-zero arguments
 */
std::vector<double> common_nonzero_args() {
  return std::vector<double>{-1.3,
                             0.49,
                             0.99,
                             1.01,
                             stan::math::positive_infinity(),
                             stan::math::negative_infinity(),
                             stan::math::not_a_number()};
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

/**
 * Test that the specified polymorphic unary function produces the
 * same results, exceptions, and has derivatives consistent with
 * finite differences as returned by the primitive version of the
 * function, when applied to the common arguments as defined by
 * `common_args()`.
 *
 * @tparam F type of polymorphic unary functor
 * @param f unary functor to test
 */
template <typename F>
void expect_common_unary(const F& f) {
  auto args = common_args();
  for (double x1 : args)
    expect_ad(f, x1);
}

/**
 * Test that the specified polymorphic unary function produces the
 * same results, exceptions, and has derivatives consistent with
 * finite differences as returned by the primitive version of the
 * function, when applied to all pairs of common arguments as defined
 * by `common_args()`.
 *
 * @tparam F type of polymorphic binary functor
 * @param f functor to test
 */
template <typename F>
void expect_common_binary(const F& f) {
  auto args = common_args();
  for (double x1 : args)
    for (double x2 : args)
      expect_ad(f, x1, x2);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has derivatives consistent with finite
 * differences as returned by the primitive version of the function
 * when applied to all common arguments as defined by `common_args`.
 *
 * <p>The function must be defined from scalars to scalars and from
 * containers to containers, always producing the same output type as
 * input type.  The value for containers must be the same as applying
 * the scalar function elementwise.
 *
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <typename F>
void expect_common_unary_vectorized(const F& f) {
  auto args = common_args();
  for (double x1 : args)
    stan::test::expect_ad_vectorized(f, x1);
}

/**
 * Test that the specified vectorized unary function produces the same
 * results and exceptions, and has derivatives consistent with finite
 * differences as returned by the primitive version of the function
 * when applied to all common non-zero arguments as defined by
 * `common_nonzero_args`.
 *
 * <p>The function must be defined from scalars to scalars and from
 * containers to containers, always producing the same output type as
 * input type.  The value for containers must be the same as applying
 * the scalar function elementwise.
 *
 * @tparam F type of functor to test
 * @param f functor to test
 */
template <typename F>
void expect_common_nonzero_unary_vectorized(const F& f) {
  auto args = common_nonzero_args();
  for (double x1 : args)
    stan::test::expect_ad_vectorized(f, x1);
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
void expect_values(const F& f, const T1& x1, const T2& x2) {
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

/**
 * Test that the specified polymorphic unary function produces the
 * same results and exceptions consistent with the primitive version
 * of the function, when applied to all pairs of common arguments as
 * defined by `common_args()`.  Derivatives are not tested, because
 * comparison operators return boolean values.
 *
 * @tparam F type of polymorphic binary functor
 * @param f functor to test
 */
template <typename F>
void expect_common_comparison(const F& f) {
  auto args = common_args();
  for (double x1 : args)
    for (double x2 : args)
      expect_values(f, x1, x2);
}

}  // namespace test
}  // namespace stan

#endif
