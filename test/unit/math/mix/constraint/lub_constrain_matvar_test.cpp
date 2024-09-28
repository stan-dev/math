#include <test/unit/math/test_ad.hpp>

namespace lub_constrain_tests {
template <typename T1, typename T2, typename T3>
void expect_matvar(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<false>(x, lb, ub, lp);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain<true>(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain<true>(x, lb, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain<true>(x, lb, ub, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad_matvar(f1, x, lb, ub);
  stan::test::expect_ad_matvar(f2, x, lb, ub);
  stan::test::expect_ad_matvar(f3, x, lb, ub);
  stan::test::expect_ad_matvar(f4, x, lb, ub);
}
}  // namespace lub_constrain_tests

TEST(mathMixMatFun, lub_constrain_scalars_matvar) {
  double x1 = 0.7;
  double x2 = -38.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);
  lub_constrain_tests::expect_matvar(x1, lb, lb);
  lub_constrain_tests::expect_matvar(x2, lb, lb);
  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub_inf);
  lub_constrain_tests::expect_matvar(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_matvar(x1, lb_inf, ub);
  lub_constrain_tests::expect_matvar(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect_matvar(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect_matvar(x2, lb_inf, ub_inf);
}

TEST(mathMixMatFun, lub_constrain_vector_scalar_scalar_matvar) {
  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);
  lub_constrain_tests::expect_matvar(x1, lb, lb);
  lub_constrain_tests::expect_matvar(x2, lb, lb);

  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub_inf);
  lub_constrain_tests::expect_matvar(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_matvar(x1, lb_inf, ub);
  lub_constrain_tests::expect_matvar(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect_matvar(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect_matvar(x2, lb_inf, ub_inf);
}

TEST(mathMixMatFun, lub_constrain_vector_vector_scalar_matvar) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 5.0, -6.0, 6.0;
  double ub = 13.5;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub_inf);
  lub_constrain_tests::expect_matvar(x2, lb, ub_inf);

  // lb inf
  lb(1, 1) = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // both inf
  lub_constrain_tests::expect_matvar(x1, lb, ub_inf);
  lub_constrain_tests::expect_matvar(x2, lb, ub_inf);
}

TEST(mathMixMatFun, lub_constrain_vector_scalar_vector_matvar) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  double lb = -2.0;
  Eigen::MatrixXd ub(2, 2);
  ub << -1.0, 5.0, 0.0, 38.0;

  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // ub inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // lb inf
  ub(1, 1) = 38.0;
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_matvar(x1, lb_inf, ub);
  lub_constrain_tests::expect_matvar(x2, lb_inf, ub);

  // both inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb_inf, ub);
  lub_constrain_tests::expect_matvar(x2, lb_inf, ub);
}

TEST(mathMixMatFun, lub_constrain_vector_vector_vector_matvar) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 3.0, -6.0, 6.0;
  Eigen::MatrixXd ub(2, 2);
  ub << -1.0, 5.0, 0.0, 38.0;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // ub inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // lb inf
  ub(1, 1) = 38.0;
  lb(1, 1) = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);

  // both inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect_matvar(x1, lb, ub);
  lub_constrain_tests::expect_matvar(x2, lb, ub);
}
