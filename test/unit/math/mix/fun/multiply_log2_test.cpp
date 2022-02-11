#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, multiplyLog2_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::multiply_log;
    return multiply_log(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;

  Eigen::VectorXd x1(3);
  x1 << 1.0, 2.0, 3.0;
  Eigen::RowVectorXd x2(3);
  x2 << 1.0, 2.0, 3.0;
  Eigen::MatrixXd x3(2, 3);
  x3 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  stan::test::expect_ad_matvar(f, x1, x1);
  stan::test::expect_ad_matvar(f, x1, 2.0);
  stan::test::expect_ad_matvar(f, 3.0, x1);
  stan::test::expect_ad_matvar(f, x2, x2);
  stan::test::expect_ad_matvar(f, x2, 2.5);
  stan::test::expect_ad_matvar(f, 3.5, x2);
  stan::test::expect_ad_matvar(f, x3, x3);
  stan::test::expect_ad_matvar(f, x3, 4.0);
  stan::test::expect_ad_matvar(f, 5.0, x3);

  Eigen::VectorXd x4(0);
  Eigen::RowVectorXd x5(0);
  Eigen::MatrixXd x6(0, 0);

  stan::test::expect_ad_matvar(f, x4, x4);
  stan::test::expect_ad_matvar(f, x4, 2.0);
  stan::test::expect_ad_matvar(f, 3.0, x4);
  stan::test::expect_ad_matvar(f, x5, x5);
  stan::test::expect_ad_matvar(f, x5, 2.5);
  stan::test::expect_ad_matvar(f, 3.5, x5);
  stan::test::expect_ad_matvar(f, x6, x6);
  stan::test::expect_ad_matvar(f, x6, 4.0);
  stan::test::expect_ad_matvar(f, 5.0, x6);
}
