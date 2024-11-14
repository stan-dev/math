#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, read_cov_matrix) {
  auto f = [](int K) {
    Eigen::VectorXd rx2 = (Eigen::VectorXd::Random(K).array() + 2.0).matrix();
    return [K, rx2](const auto& x1) {
      auto&& x1_ref = stan::math::to_ref(x1);
      stan::plain_type_t<std::decay_t<decltype(x1)>> x2
          = stan::math::add(x1_ref.head(K), rx2);
      return stan::math::read_cov_matrix(x1_ref, x2);
    };
  };

  Eigen::VectorXd x1(6);
  x1 << -0.9, 0.0, 0.99, 0.1, 0.2, 0.3;
  Eigen::VectorXd x2(3);
  x2 << -0.3, 0.2, -0.99;
  Eigen::VectorXd x3(1);
  x3 << 0.1;

  stan::test::expect_ad(f(4), x1);
  stan::test::expect_ad(f(3), x2);

  stan::test::expect_ad_matvar(f(4), x1);
  stan::test::expect_ad_matvar(f(3), x2);
}

TEST_F(mathMix, read_cov_matrix_lp) {
  auto f1 = [](int K) {
    Eigen::VectorXd rx2 = (Eigen::VectorXd::Random(K).array() + 2.0).matrix();
    return [K, rx2](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      auto&& x1_ref = stan::math::to_ref(x1);
      stan::plain_type_t<std::decay_t<decltype(x1)>> x2
          = stan::math::add(x1_ref.head(K), rx2);
      return stan::math::read_cov_matrix(x1_ref, x2, lp);
    };
  };

  auto f2 = [](int K) {
    Eigen::VectorXd rx2
        = (Eigen::VectorXd::Random(K).array() * 0.0 + 2.0).matrix();
    return [K, rx2](const auto& x1) {
      stan::scalar_type_t<decltype(x1)> lp = 0.0;
      auto&& x1_ref = stan::math::to_ref(x1);
      auto x2 = stan::math::eval(stan::math::add(x1_ref.head(K), rx2));
      stan::math::read_cov_matrix(x1_ref, x2, lp);
      return lp;
    };
  };

  Eigen::VectorXd x1(6);
  x1 << -0.9, 0.0, 0.99, 0.1, 0.2, 0.3;
  Eigen::VectorXd x2(3);
  x2 << -0.3, 0.2, -0.99;
  Eigen::VectorXd x3(1);
  x3 << 0.1;

  stan::test::expect_ad(f1(4), x1);
  stan::test::expect_ad(f1(3), x2);

  stan::test::expect_ad_matvar(f1(4), x1);
  stan::test::expect_ad_matvar(f1(3), x2);

  stan::test::expect_ad(f2(4), x1);
  stan::test::expect_ad(f2(3), x2);

  stan::test::expect_ad_matvar(f2(4), x1);
  stan::test::expect_ad_matvar(f2(3), x2);
}
