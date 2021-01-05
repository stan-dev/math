#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace lb_constrain_test {
template <typename T, typename U>
auto g1(const T& x, const U& lb) {
  return stan::math::lb_constrain(x, lb);
}
template <typename T, typename U>
auto g2(const T& x, const U& lb) {
  stan::return_type_t<T, U> lp = 0;
  return stan::math::lb_constrain(x, lb, lp);
}

template <typename T, typename U>
auto g3(const T& x, const U& lb) {
  stan::return_type_t<T, U> lp = 0;
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
  using stan::scalar_type_t;
  using stan::math::lb_constrain;
  using stan::math::promote_scalar_t;

  auto f1 = [](const auto& x, const auto& lb) {
    return lb_constrain_test::g1(x, lb);
  };

  auto f2 = [](const auto& x, const auto& lb) {
    return lb_constrain_test::g2(x, lb);
  };

  auto f3 = [](const auto& x, const auto& lb) {
    return lb_constrain_test::g3(x, lb);
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  Eigen::MatrixXd lbm(Eigen::MatrixXd::Random(2, 2));
  lbm << 1.0, -5.0, -stan::math::INFTY, 100.0;

  double lbd1 = -5.0;
  double lbd2 = -stan::math::INFTY;

  stan::test::expect_ad(f1, A, lbm);
  stan::test::expect_ad(f1, A, lbd1);
  stan::test::expect_ad(f1, A, lbd2);
  stan::test::expect_ad_matvar(f1, A, lbm);
  stan::test::expect_ad_matvar(f1, A, lbd1);
  stan::test::expect_ad_matvar(f1, A, lbd2);
  stan::test::expect_ad(f2, A, lbm);
  stan::test::expect_ad(f2, A, lbd1);
  stan::test::expect_ad(f2, A, lbd2);
  stan::test::expect_ad_matvar(f2, A, lbm);
  stan::test::expect_ad_matvar(f2, A, lbd1);
  stan::test::expect_ad_matvar(f2, A, lbd2);
  stan::test::expect_ad(f3, A, lbm);
  stan::test::expect_ad(f3, A, lbd1);
  stan::test::expect_ad(f3, A, lbd2);
  stan::test::expect_ad_matvar(f3, A, lbm);
  stan::test::expect_ad_matvar(f3, A, lbd1);
  stan::test::expect_ad_matvar(f3, A, lbd2);
}
