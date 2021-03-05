#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/offset_multiplier_constrain_helpers.hpp>

// array[] matrix, matrix, real
// array[] matrix, real, matrix
// array[] matrix, real, real
TEST(mathMixMatFun, offset_multiplier_stdvec_mat_scalar_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.1;
  Eigen::MatrixXd mu_inner(2, 3);
  mu_inner << -1.0, 1.0, -6.0, 1.0, 0.0, 0.01;
  Eigen::MatrixXd sigma_inner(2, 3);
  sigma_inner << 6.0, 3.0, 12.0, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> mu_vec{mu_inner, mu_inner};
  std::vector<Eigen::MatrixXd> sigma_vec{sigma_inner, sigma_inner};
  double mu_scal = -1.0;
  double sigma_scal = 7.0;
  offset_multiplier_constrain_tests::expect_vec(A, mu_scal, sigma_inner);
  offset_multiplier_constrain_tests::expect_vec(A, mu_inner, sigma_scal);
  offset_multiplier_constrain_tests::expect_vec(A, mu_scal, sigma_scal);
}
