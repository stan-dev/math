#ifndef TEST_UNIT_MATH_TEST_AD_HPP
#define TEST_UNIT_MATH_TEST_AD_HPP

#include <test/unit/math/is_finite.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/serializer.hpp>
#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>

namespace stan {
namespace test {

/**
 * Simple struct to hold the complete set of tolerances used to test a
 * function.  The default constructor uses default values for all
 * tolerances.  These are 1e-8 for values, 1e-4 for first derivatives,
 * 1e-3 for second derivatives, and 1e-2 for third derivatives.  The
 * names begin with the functional being evaluated, include an `fvar`
 * if the function is implemented using only forward mode, and end
 * with the quantity being calculated; for exmaple,
 * `hessian_fvar_grad_` is the gradient calculated by the `hessian`
 * function using forward-mode autodiff.
 *
 * `gradient_val_`: 1e-8;  `gradient_grad_`: 1e-4
 *
 * `gradient_fvar_val_`: 1e-8;  `gradient_fvar_grad_`: 1e-4
 *
 * `hessian_val_` : 1e-8; `hessian_grad_`: 1e-4; `hessian_hessian_`: 1e-3
 *
 * `hessian_fvar_val_` : 1e-8; `hessian_fvar_grad_`: 1e-4;
 * `hessian_fvar_hessian_`: 1e-3
 *
 * `grad_hessian_val_` : 1e-8; `grad_hessian_hessian_`: 1e-3;
 * `grad_hessian_grad_hessian_`: 1e-2
 */
struct ad_tolerances {
  double gradient_val_;
  double gradient_grad_;
  double gradient_fvar_val_;
  double gradient_fvar_grad_;
  double hesssian_val_;
  double hessian_grad_;
  double hessian_hessian_;
  double hessian_fvar_val_;
  double hessian_fvar_grad_;
  double hessian_fvar_hessian_;
  double grad_hessian_val_;
  double grad_hessian_hessian_;
  double grad_hessian_grad_hessian_;
  ad_tolerances()
      : gradient_val_(1e-8),
        gradient_grad_(1e-4),

        gradient_fvar_val_(1e-8),
        gradient_fvar_grad_(1e-4),

        hesssian_val_(1e-8),
        hessian_grad_(1e-4),
        hessian_hessian_(1e-3),

        hessian_fvar_val_(1e-8),
        hessian_fvar_grad_(1e-4),
        hessian_fvar_hessian_(1e-3),

        grad_hessian_val_(1e-8),
        grad_hessian_hessian_(1e-3),
        grad_hessian_grad_hessian_(1e-2) {}
};

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
  expect_near_rel("test_gradient fx = fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  stan::math::finite_diff_gradient_auto(f, x, fx_fd, grad_fd);
  expect_near_rel("test gradient grad_fd == grad_ad", grad_fd, grad_ad, 1e-4);
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
  expect_near_rel("gradient_fvar fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  stan::math::finite_diff_gradient_auto(f, x, fx_fd, grad_fd);
  expect_near_rel("gradient_fvar grad_fd == grad_ad", grad_fd, grad_ad, 1e-4);
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
  expect_near_rel("hessian_fvar fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  stan::math::finite_diff_hessian_auto(f, x, fx_fd, grad_fd, H_fd);
  expect_near_rel("hessian fvar grad_fd == grad_ad", grad_fd, grad_ad, 1e-4);
  expect_near_rel("hessian fvar H_fd = H_ad", H_fd, H_ad, 1e-3);
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
  expect_near_rel("hessian fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  stan::math::finite_diff_hessian_auto(f, x, fx_fd, grad_fd, H_fd);
  expect_near_rel("hessian grad_fd = grad_ad", grad_fd, grad_ad, 1e-4);
  expect_near_rel("hessian grad_fd H_fd == H_ad", H_fd, H_ad, 1e-3);
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
  expect_near_rel("grad_hessian fx == fx_ad", fx, fx_ad, 1e-8);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::MatrixXd H_fd;
  std::vector<Eigen::MatrixXd> grad_H_fd;
  stan::math::finite_diff_grad_hessian_auto(f, x, fx_fd, H_fd, grad_H_fd);
  expect_near_rel("grad hessian H_fd == H_ad", H_fd, H_ad, 1e-3);
  EXPECT_EQ(x.size(), grad_H_fd.size());
  for (size_t i = 0; i < grad_H_fd.size(); ++i)
    expect_near_rel("grad hessian grad_H_fd[i] == grad_H_ad[i]", grad_H_fd[i],
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

template <typename F, typename T>
void expect_all_throw(const F& f, const T& x) {
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
 * Tests that the function produces the same argument for the
 * specified integer and the specified integer cast to double and also
 * provides all the standard autodiff tests after casting the int to
 * double.
 *
 * Implementation note:  This is declared as an overload rather than a
 * specialization because C++1y prohibits partial template function
 * specialization.
 *
 * @tparam F type of functor to test
 * @param f functor to test
 * @param x1 first argument
 */
template <typename F>
void expect_ad_v(const F& f, const int& x1) {
  // test function produces same value on integer and cast to double
  double x1_dbl = static_cast<double>(x1);
  expect_near_relative(f(x1), f(x1_dbl));

  // test autodiff with integer cast to double
  expect_ad_v(f, x1_dbl);
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
  // neither x1 nor x2 are int
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

template <typename F, typename T1>
void expect_ad_vv(const F& f, const T1& x1, const int& x2) {
  double x2_dbl = static_cast<double>(x2);

  // test value using int arg matches value with arg cast to double
  expect_near_relative(f(x1, x2), f(x1, x2_dbl));

  // test deriv at arg cast to double
  expect_ad_vv(f, x1, x2_dbl);

  // bind int and test deriv on remainder
  auto f_x2 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x1ds = ds.read(x1);
    return serialize_return(f(x1ds, x2));
  };
  expect_ad_v(f_x2, x1);
}

template <typename F, typename T2>
void expect_ad_vv(const F& f, const int& x1, const T2& x2) {
  // test value vs. cast to double
  double x1_dbl = static_cast<double>(x1);
  expect_near_relative(f(x1, x2), f(x1_dbl, x2));

  // test deriv at arg cast to double
  expect_ad_vv(f, x1_dbl, x2);

  // bind int and test deriv on remainder
  auto f_x1 = [&](const auto& v) {
    auto ds = to_deserializer(v);
    auto x2ds = ds.read(x2);
    return serialize_return(f(x1, x2ds));
  };
  expect_ad_v(f_x1, x2);
}

template <typename F>
void expect_ad_vv(const F& f, const int& x1, const int& x2) {
  // test value vs. casts to double
  double x1_dbl = static_cast<double>(x1);
  double x2_dbl = static_cast<double>(x2);
  expect_near_relative(f(x1_dbl, x2_dbl), f(x1, x2));

  //  autodiff tests for all combinations of casts
  expect_ad_vv(f, x1, x2_dbl);
  expect_ad_vv(f, x1_dbl, x2);
  expect_ad_vv(f, x1_dbl, x2_dbl);
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

// TEST FUNCTIONS AFTER HERE
// ===================================================================

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

template <typename F>
void expect_unary_vectorized(const ad_tolerances& tols, const F& f) {}

/**
 * Base case for variadic function testing any number of arguments.
 * The base case always succeeds as there is nothing to test.
 *
 * @tparam F type of function to test
 * @param f function to test
 */
template <typename F>
void expect_unary_vectorized(const F& f) {}

template <typename F, typename T, typename... Ts>
void expect_unary_vectorized(const ad_tolerances& tols, const F& f, T x,
                             Ts... xs) {
  stan::test::expect_ad_vectorized(f, x);
  expect_unary_vectorized(tols, f, xs...);
}

/**
 * Recursive case for variadic unary vectorized function tests for the
 * specified argument and pack of arguments.
 *
 * <p>This test delegates to `expect_ad_vectorized(F, T)`; see that
 * function's documentation and requirements on the type of functions.
 *
 * @tparam F type of function to test
 * @tparam T type of first argument to test
 * @tparam Ts type of remaining arguments to test
 * @param f function to test
 * @param x argument to test
 * @param xs arguments to test
 */
template <typename F, typename T, typename... Ts>
void expect_unary_vectorized(const F& f, T x, Ts... xs) {
  ad_tolerances tols;  // default tolerances
  expect_unary_vectorized(tols, f, x, xs...);
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
void expect_comparison_vals(const F& f, const T1& x1, const T2& x2) {
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

template <typename F>
void expect_vals(const F& f, int x1) {
  using stan::math::fvar;
  using stan::math::var;
  typedef var v;
  typedef fvar<double> fd;
  typedef fvar<fvar<double>> ffd;
  typedef fvar<var> fv;
  typedef fvar<fvar<var>> ffv;

  expect_near_relative(f(x1), f(static_cast<double>(x1)));
  expect_near_relative(f(x1), f(v(x1)));
  expect_near_relative(f(x1), f(fd(x1)));
  expect_near_relative(f(x1), f(ffd(x1)));
  expect_near_relative(f(x1), f(fv(x1)));
  expect_near_relative(f(x1), f(ffv(x1)));
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
      expect_comparison_vals(f, x1, x2);
}

}  // namespace test
}  // namespace stan

#endif
