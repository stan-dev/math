#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
void expect_ub_constrain(double x, double ub) {
  auto f1 = [](const auto& x, const auto& ub) {
    return stan::math::ub_constrain(x, ub);
  };
  auto f2 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain(x, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    stan::math::ub_constrain(x, ub, lp);
    return lp;
  };
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
  using stan::scalar_type_t;
  using stan::math::promote_scalar_t;
  using stan::math::ub_constrain;

  auto f1 = [](const auto& x, const auto& ub) {
    return stan::math::ub_constrain(x, ub);
  };

  auto f2 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    return stan::math::ub_constrain(x, ub, lp);
  };

  auto f3 = [](const auto& x, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(ub)> lp = 0;
    stan::math::ub_constrain(x, ub, lp);
    return lp;
  };

  Eigen::MatrixXd A(2, 2);
  A << -1.1, 0.0, 1.0, 2.0;
  Eigen::MatrixXd ubm(2, 2);
  ubm << 1.0, -2.0, 2.0, 0.7;

  double ubd1 = -1.1;

  stan::test::expect_ad(f1, A, ubm);
  stan::test::expect_ad(f1, A, ubd1);
  stan::test::expect_ad_matvar(f1, A, ubm);
  stan::test::expect_ad_matvar(f1, A, ubd1);
  stan::test::expect_ad(f2, A, ubm);
  stan::test::expect_ad(f2, A, ubd1);
  stan::test::expect_ad_matvar(f2, A, ubm);
  stan::test::expect_ad_matvar(f2, A, ubd1);
  stan::test::expect_ad(f3, A, ubm);
  stan::test::expect_ad(f3, A, ubd1);
  stan::test::expect_ad_matvar(f3, A, ubm);
  stan::test::expect_ad_matvar(f3, A, ubd1);
}
