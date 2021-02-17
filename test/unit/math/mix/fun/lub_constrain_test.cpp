#include <test/unit/math/test_ad.hpp>

namespace lub_constrain_tests {
template <typename T1, typename T2, typename T3>
void expect(const T1& x, const T2& lb, const T3& ub) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return stan::math::lub_constrain(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain(x, lb, ub, lp);
    return lp;
  };
  auto f4 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    auto xx = stan::math::lub_constrain(x, lb, ub, lp);
    return stan::math::add(lp, stan::math::sum(xx));
  };

  stan::test::expect_ad(f1, x, lb, ub);
  stan::test::expect_ad(f2, x, lb, ub);
  stan::test::expect_ad(f3, x, lb, ub);
  stan::test::expect_ad(f4, x, lb, ub);
}

}

TEST(mathMixMatFun, lub_constrain_scalars) {

  double x1 = 0.7;
  double x2 = -38.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);
  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect(x2, lb_inf, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_scalar_scalar) {

  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);


  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect(x2, lb_inf, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_vector_scalar) {

  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  // fvar Hessian fails at 0?
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 5.0, -6.0, 6.0;
  double ub = 13.5;
  using stan::math::var;
  using stan::math::promote_scalar_t;

  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);

  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);


  // lb inf
  lb(1, 1) =  stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_scalar_vector) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return stan::math::lub_constrain(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain(x, lb, ub, lp);
    return lp;
  };

  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  double lb = -2.0;
  Eigen::MatrixXd ub(1, 1);
  ub << 3.5;
  using stan::math::var;
  using stan::math::promote_scalar_t;

  stan::test::expect_ad(f1, x1, lb, ub);
  stan::test::expect_ad(f1, x2, lb, ub);

  stan::test::expect_ad(f2, x1, lb, ub);
  stan::test::expect_ad(f2, x2, lb, ub);

  stan::test::expect_ad(f3, x1, lb, ub);
  stan::test::expect_ad(f3, x2, lb, ub);

  // ub inf
  Eigen::MatrixXd ub_inf(1, 1);
  ub_inf << stan::math::INFTY;
  stan::test::expect_ad(f1, x1, lb, ub_inf);
  stan::test::expect_ad(f1, x2, lb, ub_inf);

  stan::test::expect_ad(f2, x1, lb, ub_inf);
  stan::test::expect_ad(f2, x2, lb, ub_inf);

  stan::test::expect_ad(f3, x1, lb, ub_inf);
  stan::test::expect_ad(f3, x2, lb, ub_inf);


  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  stan::test::expect_ad(f1, x1, lb_inf, ub);
  stan::test::expect_ad(f1, x2, lb_inf, ub);

  stan::test::expect_ad(f2, x1, lb_inf, ub);
  stan::test::expect_ad(f2, x2, lb_inf, ub);

  stan::test::expect_ad(f3, x1, lb_inf, ub);
  stan::test::expect_ad(f3, x2, lb_inf, ub);

  // both inf
  stan::test::expect_ad(f1, x1, lb_inf, ub_inf);
  stan::test::expect_ad(f1, x2, lb_inf, ub_inf);

  stan::test::expect_ad(f2, x1, lb_inf, ub_inf);
  stan::test::expect_ad(f2, x2, lb_inf, ub_inf);

  stan::test::expect_ad(f3, x1, lb_inf, ub_inf);
  stan::test::expect_ad(f3, x2, lb_inf, ub_inf);

}

TEST(mathMixMatFun, lub_constrain_vector_vector_vector) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return stan::math::lub_constrain(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain(x, lb, ub, lp);
    return lp;
  };

  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  Eigen::MatrixXd lb(1, 1);
  lb << -2.0;
  Eigen::MatrixXd ub(1, 1);
  ub << 3.5;
  using stan::math::var;
  using stan::math::promote_scalar_t;

  stan::test::expect_ad(f1, x1, lb, ub);
  stan::test::expect_ad(f1, x2, lb, ub);

  stan::test::expect_ad(f2, x1, lb, ub);
  stan::test::expect_ad(f2, x2, lb, ub);

  stan::test::expect_ad(f3, x1, lb, ub);
  stan::test::expect_ad(f3, x2, lb, ub);

  // ub inf
  Eigen::MatrixXd ub_inf(1, 1);
  ub_inf << stan::math::INFTY;
  stan::test::expect_ad(f1, x1, lb, ub_inf);
  stan::test::expect_ad(f1, x2, lb, ub_inf);

  stan::test::expect_ad(f2, x1, lb, ub_inf);
  stan::test::expect_ad(f2, x2, lb, ub_inf);

  stan::test::expect_ad(f3, x1, lb, ub_inf);
  stan::test::expect_ad(f3, x2, lb, ub_inf);


  // lb inf
  Eigen::MatrixXd lb_inf(1, 1);
  lb_inf << stan::math::NEGATIVE_INFTY;
  stan::test::expect_ad(f1, x1, lb_inf, ub);
  stan::test::expect_ad(f1, x2, lb_inf, ub);

  stan::test::expect_ad(f2, x1, lb_inf, ub);
  stan::test::expect_ad(f2, x2, lb_inf, ub);

  stan::test::expect_ad(f3, x1, lb_inf, ub);
  stan::test::expect_ad(f3, x2, lb_inf, ub);

  // both inf
  stan::test::expect_ad(f1, x1, lb_inf, ub_inf);
  stan::test::expect_ad(f1, x2, lb_inf, ub_inf);

  stan::test::expect_ad(f2, x1, lb_inf, ub_inf);
  stan::test::expect_ad(f2, x2, lb_inf, ub_inf);

  stan::test::expect_ad(f3, x1, lb_inf, ub_inf);
  stan::test::expect_ad(f3, x2, lb_inf, ub_inf);

}
