#include <test/unit/math/test_ad.hpp>

// log(1 - tanh(x)^2) and its derivative -2 tanh(x) at large |x|, against
// mpmath. Forming 1 - tanh(x)^2 loses all digits from |x| = 19.06.
TEST(mathMixMatFun, corr_constrain_lp_large_argument) {
  using stan::math::fvar;
  using stan::math::var;
  Eigen::VectorXd x(3);
  x << 15.0, -19.5, 30.0;
  const double lp_ref = -124.84111691664052;
  Eigen::VectorXd grad_ref(3);
  grad_ref << -1.9999999999996257, 2.0, -2.0;

  double lp_d = 0;
  stan::math::corr_constrain(x, lp_d);
  EXPECT_NEAR(lp_ref, lp_d, 1e-13);

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_v = stan::math::to_var(x);
  var lp_v = 0;
  stan::math::corr_constrain(x_v, lp_v);
  lp_v.grad();
  EXPECT_NEAR(lp_ref, lp_v.val(), 1e-13);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad_ref(i), x_v(i).adj(), 1e-14);
  }
  stan::math::recover_memory();

  stan::math::var_value<Eigen::VectorXd> x_m(x);
  var lp_m = 0;
  stan::math::corr_constrain(x_m, lp_m);
  lp_m.grad();
  EXPECT_NEAR(lp_ref, lp_m.val(), 1e-13);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad_ref(i), x_m.adj()(i), 1e-14);
  }
  stan::math::recover_memory();

  for (int n = 0; n < 3; ++n) {
    Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> x_f(3);
    for (int i = 0; i < 3; ++i) {
      x_f(i) = fvar<double>(x(i), i == n ? 1.0 : 0.0);
    }
    fvar<double> lp_f = 0;
    stan::math::corr_constrain(x_f, lp_f);
    EXPECT_NEAR(lp_ref, lp_f.val_, 1e-13);
    EXPECT_NEAR(grad_ref(n), lp_f.d_, 1e-14);
  }
}

TEST(mathMixMatFun, corr_constrain) {
  auto f = [](const auto& x1) {
    stan::scalar_type_t<decltype(x1)> lp = 0.0;
    return stan::math::corr_constrain<false>(x1, lp);
  };

  std::vector<double> x0 = {-1.0, 2.0, 3.0};
  Eigen::VectorXd x1(3);
  x1 << -1.0, -2.0, 0.3;
  Eigen::RowVectorXd x2(3);
  x2 << -1.0, -2.0, 0.3;
  stan::test::expect_ad(f, x0);
  stan::test::expect_ad(f, x1);
  stan::test::expect_ad(f, x2);
  stan::test::expect_ad_matvar(f, x1);
  stan::test::expect_ad_matvar(f, x2);

  std::vector<double> x4;
  Eigen::VectorXd x5(0);
  Eigen::RowVectorXd x6(0);

  stan::test::expect_ad(f, x4);
  stan::test::expect_ad(f, x5);
  stan::test::expect_ad(f, x6);

  stan::test::expect_ad_matvar(f, x5);
  stan::test::expect_ad_matvar(f, x6);
}

TEST(mathMixMatFun, corr_constrain_lp) {
  auto f1 = [](const auto& x1) {
    stan::scalar_type_t<decltype(x1)> lp = 0.0;
    return stan::math::corr_constrain<true>(x1, lp);
  };

  auto f2 = [](const auto& x1) {
    stan::scalar_type_t<decltype(x1)> lp = 0.0;
    stan::math::corr_constrain<true>(x1, lp);
    return lp;
  };

  std::vector<double> x0 = {-1.0, 2.0, 3.0};
  Eigen::VectorXd x1(3);
  x1 << -1.0, -2.0, 0.3;
  Eigen::RowVectorXd x2(3);
  x2 << -1.0, -2.0, 0.3;
  stan::test::expect_ad(f1, x0);
  stan::test::expect_ad(f1, x1);
  stan::test::expect_ad(f1, x2);
  stan::test::expect_ad_matvar(f1, x1);
  stan::test::expect_ad_matvar(f1, x2);

  stan::test::expect_ad(f2, x0);
  stan::test::expect_ad(f2, x1);
  stan::test::expect_ad(f2, x2);
  stan::test::expect_ad_matvar(f2, x1);
  stan::test::expect_ad_matvar(f2, x2);

  std::vector<double> x4;
  Eigen::VectorXd x5(0);
  Eigen::RowVectorXd x6(0);

  stan::test::expect_ad(f1, x4);
  stan::test::expect_ad(f1, x5);
  stan::test::expect_ad(f1, x6);
  stan::test::expect_ad_matvar(f1, x5);
  stan::test::expect_ad_matvar(f1, x6);

  stan::test::expect_ad(f2, x4);
  stan::test::expect_ad(f2, x5);
  stan::test::expect_ad(f2, x6);
  stan::test::expect_ad_matvar(f2, x5);
  stan::test::expect_ad_matvar(f2, x6);
}
