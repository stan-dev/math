#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/constraint/offset_multiplier_constrain_helpers.hpp>

// array[] matrix, matrix, array[] matrix
// array[] matrix, matrix, matrix
TEST(mathMixMatFun, offset_multiplier_stdvec_mu_mat_sigma_vec_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.1;
  Eigen::MatrixXd mu_inner(2, 3);
  mu_inner << -1.0, 1.0, -6.0, 1.0, 0.0, 0.01;
  Eigen::MatrixXd sigma_inner(2, 3);
  sigma_inner << 6.0, 3.0, 12.0, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> mu_vec{mu_inner, mu_inner};
  std::vector<Eigen::MatrixXd> sigma_vec{sigma_inner, sigma_inner};
  offset_multiplier_constrain_tests::expect_vec(A, mu_inner, sigma_vec);
  offset_multiplier_constrain_tests::expect_vec(A, mu_inner, sigma_inner);
}

// array[] matrix, array[] matrix, array[] matrix
// array[] matrix, array[] matrix, matrix
TEST(mathMixMatFun, offset_multiplier_stdvec_mu_vec_sigma_mat_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0, 0.1;
  Eigen::MatrixXd mu_inner(2, 3);
  mu_inner << -1.0, 1.0, -6.0, 1.0, 0.0, 0.01;
  Eigen::MatrixXd sigma_inner(2, 3);
  sigma_inner << 6.0, 3.0, 12.0, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> mu_vec{mu_inner, mu_inner};
  std::vector<Eigen::MatrixXd> sigma_vec{sigma_inner, sigma_inner};
  offset_multiplier_constrain_tests::expect_vec(A, mu_vec, sigma_vec);
  offset_multiplier_constrain_tests::expect_vec(A, mu_vec, sigma_inner);
}
