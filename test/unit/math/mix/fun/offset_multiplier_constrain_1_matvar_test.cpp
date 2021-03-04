#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/offset_multiplier_constrain_matvar_helpers.hpp>

TEST(mathMixMatFun, offset_multiplier_constrain_matvar_scalars) {
  double x1 = 0.7;
  double x2 = -38.1;
  double mu = -2.0;
  double sigma = 3.5;
  offset_multiplier_constrain_tests::expect_matvar(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect_matvar(x2, mu, sigma);
}

TEST(mathMixMatFun, offset_multiplier_constrain_matvar_vector_scalar_scalar) {
  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  double mu = -2.0;
  double sigma = 3.5;
  offset_multiplier_constrain_tests::expect_matvar(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect_matvar(x2, mu, sigma);
}

TEST(mathMixMatFun, offset_multiplier_constrain_matvar_vector_vector_scalar) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd mu(2, 2);
  mu << -3.0, 5.0, -6.0, 6.0;
  double sigma = 13.5;
  offset_multiplier_constrain_tests::expect_matvar(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect_matvar(x2, mu, sigma);
}
