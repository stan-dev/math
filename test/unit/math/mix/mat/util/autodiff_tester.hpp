#ifndef TEST_UNIT_MATH_MIX_MAT_UTIL_AUTODIFF_TESTER_HPP
#define TEST_UNIT_MATH_MIX_MAT_UTIL_AUTODIFF_TESTER_HPP

#include <stan/math/mix/mat.hpp>
#include <test/unit/math/mix/mat/seq_reader.hpp>
#include <test/unit/math/mix/mat/seq_writer.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace stan {
namespace math {
namespace test {

// For an example of how to use this tester, see
// test/unit/math/mix/core/operator_addition_test.cpp

/**
 * Return true if the specified value is finite.
 *
 * @param x value to test
 * @return true if value is finite
 */
bool is_finite(double x) { return !is_inf(x) && !is_nan(x); }

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
 * Test that scalars x1 and x2 are within the specified
 * tolerance, with identity behavior for infinite and NaN
 * values.
 */
template <typename T1, typename T2>
void expect_near(const std::string& msg, const T1& x1, const T2& x2,
                 double tol = 1e-9) {
  if (is_nan(x1) || is_nan(x2))
    EXPECT_TRUE(is_nan(x1) && is_nan(x2))
        << "expect_near(" << x1 << ", " << x2 << ")" << std::endl
        << msg << std::endl;
  else if (is_inf(x1) || is_inf(x2))
    EXPECT_EQ(x1, x2) << "expect_near(" << x1 << ", " << x2 << ")" << std::endl
                      << msg << std::endl;
  else
    EXPECT_NEAR(x1, x2, tol)
        << "expect_near(" << x1 << ", " << x2 << ")" << std::endl
        << msg << std::endl;
}

/**
 * Tests that matrices (or vectors) x1 and x2 are same size and
 * have near values up to specified tolerance.
 */
template <typename T, int R, int C>
void expect_near(const std::string& msg, const Eigen::Matrix<T, R, C>& x1,
                 const Eigen::Matrix<T, R, C>& x2, double tol = 1e-7) {
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
 * Tests that the function f applied to the argument x yields
 * the expected value fx (with scalars as double).
 */
template <typename F>
void test_value(const F& f, const Eigen::VectorXd& x, double fx) {
  if (is_nan(fx))
    EXPECT_TRUE(is_nan(f(x))) << "test_value is_nan(" << f(x) << std::endl;
  else
    expect_near("test_value fx == f(x)", fx, f(x));
}

/**
 * Tests that the function f applied to the argument x yields
 * the value fx and the correct first-order derivatives as
 * calaculated with the gradient functional using var.
 */
template <typename F>
void test_gradient(const F& f, const Eigen::VectorXd& x, double fx,
                   bool test_derivs) {
  Eigen::VectorXd grad_ad;
  double fx_ad = fx;
  gradient<F>(f, x, fx_ad, grad_ad);
  expect_near("test_gradient fx = fx_ad", fx, fx_ad);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  finite_diff_gradient(f, x, fx_fd, grad_fd);
  expect_near("test gradient grad_fd == grad_ad", grad_fd, grad_ad);
}

/**
 * Tests that the function f applied to the argument x yields
 * the value fx and correct first-order derivatives as
 * calculated by the gradient functionional using fvar<double>
 * scalars.
 */
template <typename F>
void test_gradient_fvar(const F& f, const Eigen::VectorXd& x, double fx,
                        bool test_derivs) {
  Eigen::VectorXd grad_ad;
  double fx_ad = fx;
  gradient<double, F>(f, x, fx_ad, grad_ad);
  expect_near("gradient_fvar fx == fx_ad", fx, fx_ad);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  Eigen::VectorXd grad_fd;
  double fx_fd;
  finite_diff_gradient(f, x, fx_fd, grad_fd);
  expect_near("gradeint_fvar gard_fd == grad_ad", grad_fd, grad_ad);
}

/**
 * Tests that the function f applied to the argument x yields
 * the value fx and correct first- and second-order derivatives
 * as calculated by the hessian functional using fvar<var>
 * scalars.
 */
template <typename F>
void test_hessian_fvar(const F& f, const Eigen::VectorXd& x, double fx,
                       bool test_derivs) {
  double fx_ad;
  Eigen::VectorXd grad_ad;
  Eigen::MatrixXd H_ad;
  hessian<double, F>(f, x, fx_ad, grad_ad, H_ad);
  expect_near("hessian_fvar fx == fx_ad", fx, fx_ad);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
  expect_near("hessian fvar grad_fd == grad_ad", grad_fd, grad_ad);
  expect_near("hessian fvar H_fd = H_ad", H_fd, H_ad);
}

/**
 * Tests that the function f applied to the argument x yields
 * the value fx and correct first- and second-order derivatives
 * as calculated by the hessian functional using
 * fvar<fvar<double>> scalars.
 */
template <typename F>
void test_hessian(const F& f, const Eigen::VectorXd& x, double fx,
                  bool test_derivs) {
  double fx_ad;
  Eigen::VectorXd grad_ad;
  Eigen::MatrixXd H_ad;
  hessian<F>(f, x, fx_ad, grad_ad, H_ad);
  expect_near("hessian fx == fx_ad", fx, fx_ad);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::VectorXd grad_fd;
  Eigen::MatrixXd H_fd;
  finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
  expect_near("hessian grad_fd = grad_ad", grad_fd, grad_ad);
  expect_near("hessian grad_fd H_fd == H_ad", H_fd, H_ad);
}

/**
 * Tests that the function f applied to the argument x yields
 * the value fx and correct first-, second-, and third-order
 * derivatives as calculated by the hessian functional using
 * fvar<fvar<var>> scalars.
 */
template <typename F>
void test_grad_hessian(const F& f, const Eigen::VectorXd& x, double fx,
                       bool test_derivs) {
  double fx_ad;
  Eigen::MatrixXd H_ad;
  std::vector<Eigen::MatrixXd> grad_H_ad;
  grad_hessian(f, x, fx_ad, H_ad, grad_H_ad);
  expect_near("grad_hessian fx == fx_ad", fx, fx_ad);
  if (!test_derivs || !is_finite(x) || !is_finite(fx))
    return;
  double fx_fd;
  Eigen::MatrixXd H_fd;
  std::vector<Eigen::MatrixXd> grad_H_fd;
  finite_diff_grad_hessian(f, x, fx_fd, H_fd, grad_H_fd);
  expect_near("grad hessian H_fd == H_ad", H_fd, H_ad);
  EXPECT_EQ(x.size(), grad_H_fd.size());
  for (size_t i = 0; i < grad_H_fd.size(); ++i)
    expect_near("grad hessian grad_H_fd[i] == grad_H_ad[i]", grad_H_fd[i],
                grad_H_ad[i]);
}

template <typename T, typename F>
void expect_throw(const F& f, const Eigen::VectorXd& x) {
  Eigen::Matrix<T, -1, 1> x_t(x.rows());
  for (int i = 0; i < x.rows(); ++i)
    x_t(i) = x(i);
  try {
    f(x_t);
    FAIL() << "double throws, expect autodiff version to throw";
  } catch (...) {
    SUCCEED();
  }
}

template <typename F>
void expect_all_throw(const F& f, const Eigen::VectorXd& x) {
  using stan::math::fvar;
  using stan::math::var;
  expect_throw<var>(f, x);
  expect_throw<fvar<double> >(f, x);
  expect_throw<fvar<fvar<double> > >(f, x);
  expect_throw<fvar<var> >(f, x);
  expect_throw<fvar<fvar<var> > >(f, x);
}

// test value and derivative in all functionals vs. finite diffs
template <typename F>
void test_functor(const F& f, const Eigen::VectorXd& x, double fx,
                  bool test_derivs, bool expect_exception) {
  if (expect_exception) {
    expect_all_throw(f, x);
    return;
  }
  test_value(f, x, fx);
  test_gradient(f, x, fx, test_derivs);
  test_gradient_fvar(f, x, fx, test_derivs);
  test_hessian(f, x, fx, test_derivs);
  test_hessian_fvar(f, x, fx, test_derivs);
  test_grad_hessian(f, x, fx, test_derivs);
}

/**
 * Structure to adapt a two-argument function specified as a
 * class with a static apply(, ) method that returns a scalar to
 * a function that operates on an Eigen vector with templated
 * scalar type and returns a scalar of the same type.
 *
 * <p>It works by adapting the two-argument function to be a
 * vector function with zero, one or both arguments being
 * instantiated to doubles and the remaining arguments being
 * templated on the functor argument scalar type.
 *
 * @tparam F class with static apply(, ) method
 * @tparam T1 type of first argument with double scalars
 * @tparam T2 type of second argument with double scalars
 */
template <typename F, typename T1, typename T2>
struct binder_binary {
  T1 x1_;
  T2 x2_;
  bool fixed1_;
  bool fixed2_;

  binder_binary(const T1& x1, const T2& x2)
      : x1_(x1), x2_(x2), fixed1_(false), fixed2_(false) {}

  template <typename T>
  T operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    if (fixed1_ && fixed2_) {
      return F::apply(x1_, x2_);
    } else if (fixed1_ && !fixed2_) {
      seq_reader<T> r(theta);
      return F::apply(x1_, r.read(x2_));
    } else if (!fixed1_ && fixed2_) {
      seq_reader<T> r(theta);
      return F::apply(r.read(x1_), x2_);
    } else if (!fixed1_ && !fixed2_) {
      seq_reader<T> r(theta);
      T read_x1 = r.read(x1_);
      T read_x2 = r.read(x2_);
      return F::apply(read_x1, read_x2);
    }
    throw std::logic_error("binder_binary illegal state");
  }
};

/**
 * Tests whether the binary function specified through the static
 * function F::apply with arguments of the specified types
 * instantiated with double scalars returns right values and
 * first, second, and third derivatives.  It tests all possible
 * combinations of double, var, fvar<double>,
 * fvar<fvar<double>>, fvar<var>, and fvar<fvar<var>>
 * instantiations.
 */
template <typename F, typename T1, typename T2>
void test_ad(const T1& x1, const T2& x2, double fx, bool test_derivs,
             bool expect_exception) {
  // create binder then test all autodiff/double combos
  binder_binary<F, T1, T2> f(x1, x2);

  // test (double, double) instantiation
  f.fixed1_ = true;
  f.fixed2_ = true;
  seq_writer<double> a;
  test_functor(f, a.vector(), fx, test_derivs, expect_exception);

  // test (double, autodiff) instantiation
  f.fixed1_ = true;
  f.fixed2_ = false;
  seq_writer<double> b;
  b.write(x2);
  test_functor(f, b.vector(), fx, test_derivs, expect_exception);

  // test (autodiff, double) instantiation
  f.fixed1_ = false;
  f.fixed2_ = true;
  seq_writer<double> c;
  c.write(x1);
  test_functor(f, c.vector(), fx, test_derivs, expect_exception);

  // test (autodiff, autodiff) instantiation
  f.fixed1_ = false;
  f.fixed2_ = false;
  seq_writer<double> d;
  d.write(x1);
  d.write(x2);
  test_functor(f, d.vector(), fx, test_derivs, expect_exception);
}

template <typename F, bool is_comparison>
void test_args(const std::vector<double>& xs1, const std::vector<double>& xs2) {
  using stan::math::test::test_ad;

  // avoid testing derivatives for comparisons of equal values
  // as results will be non-differentiable in one direction
  for (size_t i = 0; i < xs1.size(); ++i) {
    for (size_t j = 0; j < xs2.size(); ++j) {
      double fx;
      bool threw_exception = false;
      try {
        fx = F::apply(xs1[i], xs2[j]);
      } catch (...) {
        threw_exception = true;
      }
      test_ad<F>(xs1[i], xs2[j], fx, !(is_comparison && xs1[i] == xs2[j]),
                 threw_exception);
    }
  }
}

template <typename F, bool is_comparison>
void test_args(const std::vector<double>& xs) {
  test_args<F, is_comparison>(xs, xs);
}

template <typename F, bool is_comparison>
void test_common_args() {
  std::vector<double> xs;
  xs.push_back(0.5);
  xs.push_back(0);
  xs.push_back(-1.3);
  xs.push_back(stan::math::positive_infinity());
  xs.push_back(stan::math::negative_infinity());
  xs.push_back(stan::math::not_a_number());
  test_args<F, is_comparison>(xs, xs);
}

}  // namespace test
}  // namespace math
}  // namespace stan
#endif
