#ifndef STAN_TEST_UNIT_MATH_MIX_CONSTRAINT_CORR_TRANSFORM_HELPERS_HPP
#define STAN_TEST_UNIT_MATH_MIX_CONSTRAINT_CORR_TRANSFORM_HELPERS_HPP

#include <test/unit/math/test_ad.hpp>
#include <cmath>

namespace corr_transform_test {

/**
 * Return the lower triangle of a matrix, row by row, without (`strict`) or
 * with the diagonal.
 */
template <typename M>
Eigen::Matrix<stan::value_type_t<M>, Eigen::Dynamic, 1> lower_triangle(
    const M& m, bool strict) {
  const Eigen::Index K = m.rows();
  Eigen::Matrix<stan::value_type_t<M>, Eigen::Dynamic, 1> v(
      strict ? (K * (K - 1)) / 2 : (K * (K + 1)) / 2);
  Eigen::Index k = 0;
  for (Eigen::Index i = 0; i < K; ++i) {
    for (Eigen::Index j = 0; j < (strict ? i : i + 1); ++j) {
      v(k++) = m(i, j);
    }
  }
  return v;
}

/**
 * Functor for `stan::math::jacobian`: the lower triangle of `constrain(y)`.
 */
template <typename F>
struct lower_triangle_of {
  F constrain;
  bool strict;
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& y) const {
    return lower_triangle(constrain(y), strict);
  }
};

/**
 * Expect `lp` from `constrain_lp(y, lp)` to be the log absolute determinant
 * of the Jacobian of the map from `y` to the lower triangle of
 * `constrain(y)`.
 */
template <typename F, typename G>
void expect_lp_is_log_det_jacobian(F constrain, G constrain_lp,
                                   const Eigen::VectorXd& y, bool strict,
                                   double tol) {
  Eigen::VectorXd fx;
  Eigen::MatrixXd J;
  stan::math::jacobian(lower_triangle_of<F>{constrain, strict}, y, fx, J);
  double lp = 0;
  constrain_lp(y, lp);
  EXPECT_NEAR(std::log(std::fabs(J.determinant())), lp, tol);
}

/**
 * Expect `lp` from `constrain_lp(y, lp)` and its gradient to match the
 * references in double, `Eigen::Matrix<var, -1, 1>`,
 * `var_value<Eigen::VectorXd>` and `fvar<double>` (one direction at a time).
 */
template <typename G>
void expect_lp_all_modes(G constrain_lp, const Eigen::VectorXd& y,
                         double lp_ref, const Eigen::VectorXd& grad_ref,
                         double tol_lp, double tol_grad) {
  using stan::math::fvar;
  using stan::math::var;
  const Eigen::Index n = y.size();

  double lp_d = 0;
  constrain_lp(y, lp_d);
  EXPECT_NEAR(lp_ref, lp_d, tol_lp) << "double";

  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  var lp_v = 0;
  constrain_lp(y_v, lp_v);
  lp_v.grad();
  EXPECT_NEAR(lp_ref, lp_v.val(), tol_lp) << "var";
  for (Eigen::Index i = 0; i < n; ++i) {
    EXPECT_NEAR(grad_ref(i), y_v(i).adj(), tol_grad) << "var, element " << i;
  }
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> y_m(y);
  var lp_m = 0;
  constrain_lp(y_m, lp_m);
  lp_m.grad();
  EXPECT_NEAR(lp_ref, lp_m.val(), tol_lp) << "var_value";
  for (Eigen::Index i = 0; i < n; ++i) {
    EXPECT_NEAR(grad_ref(i), y_m.adj()(i), tol_grad)
        << "var_value, element " << i;
  }
  stan::math::recover_memory();

  for (Eigen::Index d = 0; d < n; ++d) {
    Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> y_f(n);
    for (Eigen::Index i = 0; i < n; ++i) {
      y_f(i) = fvar<double>(y(i), i == d ? 1.0 : 0.0);
    }
    fvar<double> lp_f = 0;
    constrain_lp(y_f, lp_f);
    EXPECT_NEAR(lp_ref, lp_f.val_, tol_lp) << "fvar, direction " << d;
    EXPECT_NEAR(grad_ref(d), lp_f.d_, tol_grad) << "fvar, direction " << d;
  }
}

/**
 * Expect the lower triangle of `constrain(y)` to match `ref` within
 * `abs_tol + rel_tol * |ref|` in double, `Eigen::Matrix<var, -1, 1>`,
 * `var_value<Eigen::VectorXd>` and `fvar<double>`.
 */
template <typename F>
void expect_values_all_modes(F constrain, const Eigen::VectorXd& y, bool strict,
                             const Eigen::VectorXd& ref, double abs_tol,
                             double rel_tol) {
  using stan::math::fvar;
  using stan::math::var;
  auto expect_ref = [&](const Eigen::VectorXd& v, const char* mode) {
    ASSERT_EQ(ref.size(), v.size()) << mode;
    for (Eigen::Index i = 0; i < ref.size(); ++i) {
      EXPECT_NEAR(ref(i), v(i), abs_tol + rel_tol * std::fabs(ref(i)))
          << mode << ", element " << i;
    }
  };
  expect_ref(lower_triangle(stan::math::value_of(constrain(y)), strict),
             "double");
  Eigen::Matrix<var, Eigen::Dynamic, 1> y_v = stan::math::to_var(y);
  expect_ref(lower_triangle(stan::math::value_of(constrain(y_v)), strict),
             "var");
  stan::math::var_value<Eigen::VectorXd> y_m(y);
  expect_ref(lower_triangle(stan::math::value_of(constrain(y_m)), strict),
             "var_value");
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> y_f(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    y_f(i) = fvar<double>(y(i), 1.0);
  }
  expect_ref(lower_triangle(stan::math::value_of(constrain(y_f)), strict),
             "fvar");
  stan::math::recover_memory();
}

}  // namespace corr_transform_test

#endif
