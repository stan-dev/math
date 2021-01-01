#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
template <typename T, typename U>
stan::return_type_t<T, U> g1(const T& x, const U& ub) {
  return stan::math::ub_constrain(x, ub);
}
template <typename T, typename U>
stan::return_type_t<T, U> g2(const T& x, const U& ub) {
  T lp = 0;
  return stan::math::ub_constrain(x, ub, lp);
}
template <typename T, typename U>
stan::return_type_t<T, U> g3(const T& x, const U& ub) {
  T lp = 0;
  stan::math::ub_constrain(x, ub, lp);
  return lp;
}

void expect_ub_constrain(double x, double ub) {
  auto f1 = [](const auto& x, const auto& ub) { return g1(x, ub); };
  auto f2 = [](const auto& x, const auto& ub) { return g2(x, ub); };
  auto f3 = [](const auto& x, const auto& ub) { return g3(x, ub); };
  stan::test::expect_ad(f1, x, ub);
  stan::test::expect_ad(f2, x, ub);
  stan::test::expect_ad(f3, x, ub);
}
}  // namespace ub_constrain_test

TEST(mathMixScalFun, ub_constrain) {
  ub_constrain_test::expect_ub_constrain(-1, 2);
  ub_constrain_test::expect_ub_constrain(2, 4);
}

TEST(mathMixMatFun, ub_mat_constrain) {
  using stan::math::promote_scalar_t;
  using stan::scalar_type_t;
  using stan::math::ub_constrain;
  auto tester = [](const auto& x) {
    return ub_constrain(x, 3.0);
  };

  auto promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    x_scalar ub = 3.0;
    return ub_constrain(x, ub);
  };

  auto multi_promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    using ub_mat_t = Eigen::Matrix<x_scalar, -1, -1>;
    x_scalar ub = 3.0;
    ub_mat_t ub_mat = ub_mat_t::Constant(x.rows(), x.cols(), ub);
    return ub_constrain(x, ub_mat);
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  stan::test::expect_ad(tester, A);
  stan::test::expect_ad_matvar(tester, A);
  stan::test::expect_ad(promote_tester, A);
  stan::test::expect_ad_matvar(promote_tester, A);
  stan::test::expect_ad_matvar(multi_promote_tester, A);
}

TEST(mathMixMatFun, ub_lp_mat_constrain) {
  using stan::math::promote_scalar_t;
  using stan::scalar_type_t;
  using stan::math::ub_constrain;
  auto ub_lp_promote_tester = [](const auto& x) {
    using scalar_x = scalar_type_t<decltype(x)>;
    scalar_x ub = 3.0;
    scalar_x lp = 0;
    stan::math::ub_constrain(x, ub, lp);
    return lp;
  };
  auto multi_ub_lp_promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    using ub_mat_t = Eigen::Matrix<x_scalar, -1, -1>;
    x_scalar ub = 3.0;
    ub_mat_t ub_mat = ub_mat_t::Constant(x.rows(), x.cols(), ub);
    x_scalar lp = 0;
    return ub_constrain(x, ub_mat, lp);
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  stan::test::expect_ad(ub_lp_promote_tester, A);
  stan::test::expect_ad_matvar(ub_lp_promote_tester, A);

  stan::test::expect_ad(multi_ub_lp_promote_tester, A);
  stan::test::expect_ad_matvar(multi_ub_lp_promote_tester, A);
}
