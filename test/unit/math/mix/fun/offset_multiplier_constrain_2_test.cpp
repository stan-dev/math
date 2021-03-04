#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/offset_multiplier_constrain_helpers.hpp>

TEST(mathMixMatFun, offset_multiplier_constrain_vector_scalar_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  double mu = -2.0;
  Eigen::MatrixXd sigma(2, 2);
  sigma << -1.0, 5.0, 0.0, 38.0;

  offset_multiplier_constrain_tests::expect(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect(x2, mu, sigma);
}

TEST(mathMixMatFun, offset_multiplier_constrain_vector_vector_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.0000001, 1.0, 3.0;
  Eigen::MatrixXd mu(2, 2);
  mu << -3.0, 0.0, -6.0, 6.0;
  Eigen::MatrixXd sigma(2, 2);
  sigma << -1.0, 5.0, 0.0, 38.0;
  offset_multiplier_constrain_tests::expect(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect(x2, mu, sigma);
}
