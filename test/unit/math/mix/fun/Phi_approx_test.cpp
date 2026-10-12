#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

TEST(mathMixMatFun, PhiApprox) {
  auto f = [](const auto& x1) { return stan::math::Phi_approx(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -3.0, 1, 1.3, 3);
}

TEST(mathMixMatFun, PhiApprox_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::Phi_approx;
    return Phi_approx(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-3.0, 1, 1.3, 3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

namespace Phi_approx_test {
// d/dx Phi_approx(x) with PHI_APPROX_CUBIC_COEF and PHI_APPROX_LINEAR_COEF, by
// mpmath at 60 digits; at x = 40
// and -40 the true value, about 4e-1987, underflows to 0
inline std::vector<double> args() {
  return {-40, -10, -8, -7, -1, 0, 1, 7, 8, 10, 40};
}
inline std::vector<double> derivs() {
  return {0.0,
          5.934767576595731e-37,
          8.691773531236095e-21,
          5.126988948527854e-15,
          0.2415281648469019,
          0.399400718560196,
          0.2415281648469019,
          5.126988948527854e-15,
          8.691773531236095e-21,
          5.934767576595731e-37,
          0.0};
}
// Relative check, absolute below the smallest normal double.
inline void expect_rel(double expected, double actual, double x,
                       double tol = 1e-13) {
  EXPECT_NEAR(
      expected, actual,
      tol * std::fmax(std::fabs(expected), std::numeric_limits<double>::min()))
      << "x = " << x;
}
}  // namespace Phi_approx_test

TEST(mathMixMatFun, PhiApproxScalarDerivativeTails) {
  using stan::math::var;
  std::vector<double> xs = Phi_approx_test::args();
  std::vector<double> ds = Phi_approx_test::derivs();
  for (size_t i = 0; i < xs.size(); ++i) {
    var a = xs[i];
    stan::math::Phi_approx(a).grad();
    Phi_approx_test::expect_rel(ds[i], a.adj(), xs[i]);
    stan::math::recover_memory();
  }
}

TEST(mathMixMatFun, PhiApproxDerivativeTails) {
  using stan::math::Phi_approx;
  using stan::math::sum;
  using stan::math::var;
  std::vector<double> xs = Phi_approx_test::args();
  std::vector<double> ds = Phi_approx_test::derivs();
  Eigen::VectorXd x = stan::math::to_vector(xs);

  stan::math::var_value<Eigen::VectorXd> xv(x);
  sum(Phi_approx(xv)).grad();
  for (int i = 0; i < x.size(); ++i) {
    Phi_approx_test::expect_rel(ds[i], xv.adj()(i), xs[i]);
  }
  stan::math::recover_memory();

  Eigen::Matrix<var, Eigen::Dynamic, 1> xm = x;
  sum(Phi_approx(xm)).grad();
  for (int i = 0; i < x.size(); ++i) {
    Phi_approx_test::expect_rel(ds[i], xm(i).adj(), xs[i]);
  }
  stan::math::recover_memory();

  std::vector<var> xs_v(xs.begin(), xs.end());
  sum(Phi_approx(xs_v)).grad();
  for (size_t i = 0; i < xs.size(); ++i) {
    Phi_approx_test::expect_rel(ds[i], xs_v[i].adj(), xs[i]);
  }
  stan::math::recover_memory();
}
