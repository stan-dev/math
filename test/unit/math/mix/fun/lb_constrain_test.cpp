#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace lb_constrain_test {
template <typename T, typename U>
auto g1(const T& x, const U& lb) {
  return stan::math::lb_constrain(x, lb);
}
template <typename T, typename U>
stan::return_type_t<T, U> g2(const T& x, const U& lb) {
  T lp = 0;
  return stan::math::lb_constrain(x, lb, lp);
}

template <typename T, typename U>
stan::return_type_t<T, U> g3(const T& x, const U& lb) {
  T lp = 0;
  stan::math::lb_constrain(x, lb, lp);
  return lp;
}
}  // namespace lb_constrain_test

void expect_lb_constrain(double x, double lb) {
  auto f1 = [](const auto& x, const auto& lb) {
    return lb_constrain_test::g1(x, lb);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    return lb_constrain_test::g2(x, lb);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    return lb_constrain_test::g3(x, lb);
  };
  stan::test::expect_ad(f1, x, lb);
  stan::test::expect_ad(f2, x, lb);
  stan::test::expect_ad(f3, x, lb);
}

TEST(mathMixScalFun, lbConstrain) {
  expect_lb_constrain(-1, 2);
  expect_lb_constrain(2, 4);
}


TEST(mathMixMatFun, lb_mat_constrain) {
  using stan::math::promote_scalar_t;
  using stan::scalar_type_t;
  using stan::math::lb_constrain;
  auto tester = [](const auto& x) {
    return lb_constrain(x, 3.0);
  };

  auto promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    x_scalar lb = 3.0;
    return lb_constrain(x, lb);
  };

  auto multi_promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    using lb_mat_t = Eigen::Matrix<x_scalar, -1, -1>;
    x_scalar lb = 3.0;
    lb_mat_t lb_mat = lb_mat_t::Constant(x.rows(), x.cols(), lb);
    return lb_constrain(x, lb_mat);
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  stan::test::expect_ad(tester, A);
  stan::test::expect_ad_matvar(tester, A);
  stan::test::expect_ad(promote_tester, A);
  stan::test::expect_ad_matvar(promote_tester, A);
  stan::test::expect_ad_matvar(multi_promote_tester, A);
}

TEST(mathMixMatFun, lb_lp_mat_constrain) {
  using stan::math::promote_scalar_t;
  using stan::scalar_type_t;
  using stan::math::lb_constrain;
  auto lb_lp_promote_tester = [](const auto& x) {
    using scalar_x = scalar_type_t<decltype(x)>;
    scalar_x lb = 3.0;
    scalar_x lp = 0;
    stan::math::lb_constrain(x, lb, lp);
    return lp;
  };
  auto multi_lb_lp_promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    using lb_mat_t = Eigen::Matrix<x_scalar, -1, -1>;
    x_scalar lb = 3.0;
    lb_mat_t lb_mat = lb_mat_t::Constant(x.rows(), x.cols(), lb);
    x_scalar lp = 0;
    return lb_constrain(x, lb_mat, lp);
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  stan::test::expect_ad(lb_lp_promote_tester, A);
  stan::test::expect_ad_matvar(lb_lp_promote_tester, A);

  stan::test::expect_ad(multi_lb_lp_promote_tester, A);
  stan::test::expect_ad_matvar(multi_lb_lp_promote_tester, A);
}
