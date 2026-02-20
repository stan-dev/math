#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <vector>

struct finite_diff_typed_args_functor {
  template <typename T_theta>
  auto operator()(const std::vector<T_theta>& theta,
                  const std::vector<int>& x_i) const {
    return theta[0] * x_i[0] + theta[1] * x_i[1];
  }
};

struct finite_diff_interleaved_args_functor {
  template <typename T_a, typename T_b>
  auto operator()(const std::vector<T_a>& a, const std::vector<int>& x_i,
                  const std::vector<T_b>& b, const double c) const {
    return a[0] * x_i[0] + b[0] * x_i[1] + c;
  }
};

TEST(FwdFunctor, finiteDiffTypedIntArgs) {
  using stan::math::fvar;
  std::vector<fvar<double>> theta(2);
  theta[0] = fvar<double>(1.2, 0.5);
  theta[1] = fvar<double>(-0.7, -1.0);
  const std::vector<int> x_i{3, 4};

  auto out = stan::math::finite_diff(finite_diff_typed_args_functor{}, theta,
                                     x_i);

  EXPECT_NEAR(0.8, out.val_, 1e-10);
  EXPECT_NEAR(-2.5, out.d_, 1e-8);
}

TEST(FwdFunctor, finiteDiffTypedIntArgsFvarVar) {
  using stan::math::fvar;
  using stan::math::var;
  std::vector<fvar<var>> theta(2);
  theta[0] = fvar<var>(1.2, 0.5);
  theta[1] = fvar<var>(-0.7, -1.0);
  const std::vector<int> x_i{3, 4};

  auto out = stan::math::finite_diff(finite_diff_typed_args_functor{}, theta,
                                     x_i);

  EXPECT_NEAR(0.8, out.val_.val(), 1e-10);
  EXPECT_NEAR(-2.5, out.d_.val(), 1e-8);
}

TEST(FwdFunctor, finiteDiffInterleavedTopLevelArgs) {
  using stan::math::fvar;
  std::vector<fvar<double>> a(1);
  a[0] = fvar<double>(1.2, 0.7);
  const std::vector<int> x_i{3, 4};
  std::vector<fvar<double>> b(1);
  b[0] = fvar<double>(-0.5, -1.3);
  const double c = 1.1;

  auto out
      = stan::math::finite_diff(finite_diff_interleaved_args_functor{}, a, x_i,
                                b, c);

  EXPECT_NEAR(2.7, out.val_, 1e-10);
  EXPECT_NEAR(-3.1, out.d_, 1e-8);
}

/**
 * Regression for https://github.com/stan-dev/math/issues/3280
 *
 * Before the tuple-native finite-diff refactor, the forward-mode finite_diff
 * path could serialize and deserialize arguments in a way that converted
 * `x_i` from `std::vector<int>` into floating-point vectors. That then failed
 * to match integrate_1d's required callback signature.
 *
 * This test verifies `finite_diff` + `integrate_1d` now works with mixed
 * autodiff/real/int argument packs and produces the expected derivative.
 */
struct issue_3280_integrate_1d_functor {
  template <typename T_x, typename T_xc, typename T_theta>
  auto operator()(const T_x& x, const T_xc& /* xc */,
                  const std::vector<T_theta>& theta,
                  const std::vector<double>& x_r,
                  const std::vector<int>& x_i,
                  std::ostream* /* msgs */) const {
    return theta[0] * x + x_r[0] + x_i[0];
  }
};

TEST(FwdFunctor, finiteDiffIntegrate1dIssue3280Regression) {
  using stan::math::fvar;
  const double relative_tolerance = 1e-8;
  std::ostream* msgs = nullptr;

  std::vector<fvar<double>> theta(1);
  theta[0] = fvar<double>(2.0, 1.0);
  const std::vector<double> x_r{1.5};
  const std::vector<int> x_i{2};

  auto out = stan::math::integrate_1d(issue_3280_integrate_1d_functor{}, 0.0,
                                      1.0, theta, x_r, x_i, msgs,
                                      relative_tolerance);

  EXPECT_NEAR(4.5, out.val_, 1e-8);
  EXPECT_NEAR(0.5, out.d_, 1e-6);
}
