#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, owensT) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::owens_t(x1, x2);
  };
  stan::test::expect_ad(f, -10.9, 13);
  stan::test::expect_ad(f, 0.5, 1.0);
  stan::test::expect_ad(f, 1.0, 2.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}

TEST(mathMixScalFun, owensT_varmat) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::owens_t(x1, x2);
  };
  double scal = 2.0;
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(2, 2);
  Eigen::VectorXd vec = Eigen::VectorXd::Random(2);
  stan::test::expect_ad_matvar(f, mat, mat);
  stan::test::expect_ad_matvar(f, mat, scal);
  stan::test::expect_ad_matvar(f, scal, mat);
  stan::test::expect_ad_matvar(f, vec, vec);
  stan::test::expect_ad_matvar(f, vec, scal);
  stan::test::expect_ad_matvar(f, scal, vec);
}

TEST(mathMixScalFun, owensT_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::owens_t;
    return owens_t(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 2.0, 2.0;
  Eigen::VectorXd in2(2);
  in2 << 3.0, 4.0;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

TEST(mathMixScalFun, owensT_vec_matvar) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::owens_t;
    return owens_t(x1, x2);
  };

  Eigen::MatrixXd in1(2, 2);
  in1 << 0.5, 3.4, 5.2, 0.5;
  Eigen::MatrixXd in2(2, 2);
  in2 << 3.3, 0.9, 6.7, 3.3;
  stan::test::expect_ad_vectorized_matvar(f, in1, in2);
}
