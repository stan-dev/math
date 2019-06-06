#ifndef TEST_UNIT_MATH_TEST_AD_HPP
#define TEST_UNIT_MATH_TEST_AD_HPP

#include <stan/math.hpp>
#include <test/unit/math/util.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace stan {
namespace test {

/**
 * For the specified functor and argument, test that automatic
 * differentiation provides the value as the double-based version and
 * the same derivatives as finite differences over the double version.
 * The functor must be a map from Eigen vectors to scalars of the
 * same scalar type.
 *
 * @tparam G type of polymorphic functor
 * @param g polymorphic functor from vectors to scalars
 * @param x argument to test
 */
template <typename G>
void expect_ad_derivatives(const G& g, const Eigen::VectorXd& x) {
  double gx = g(x);
  stan::math::test::test_gradient(g, x, gx);
  stan::math::test::test_gradient_fvar(g, x, gx);
  stan::math::test::test_hessian(g, x, gx);
  stan::math::test::test_hessian_fvar(g, x, gx);
  stan::math::test::test_grad_hessian(g, x, gx);
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
 * @param h serialized functor returning a single component of
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
    stan::math::test::expect_all_throw(h(0), x);
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
 * <p>Tests all three possible instantiations of autodiff variables:
 * first argument only, second argument only, and both arguments.
 * Tests autodiff levels `rev`, `fvar<double>`, `fvar<fvar<double>>`,
 * `fvar<rev>`, and `fvar<fvar<rev>>`.
 *
 * <p>Invokes Google test framework to raise error if test fails.
 *
 * @tparam F type of polymorphic functor to test
 * @tparam T1 type of double- or int-based first argument
 * @tparam T2 type of double- or int-based second argument
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

}  // namespace test
}  // namespace stan

#endif
