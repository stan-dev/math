#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/offset_multiplier_constrain_matvar_helpers.hpp>

TEST(mathMixMatFun, offset_multiplier_constrain_matvar_vector_scalar_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  double mu = -2.0;
  Eigen::MatrixXd sigma(2, 2);
  sigma << -1.0, 5.0, 0.0, 38.0;

  offset_multiplier_constrain_tests::expect_matvar(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect_matvar(x2, mu, sigma);
}

TEST(mathMixMatFun, offset_multiplier_constrain_matvar_vector_vector_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.0000001, 1.0, 3.0;
  Eigen::MatrixXd mu(2, 2);
  mu << -3.0, 0.0, -6.0, 6.0;
  Eigen::MatrixXd sigma(2, 2);
  sigma << -1.0, 5.0, 0.0, 38.0;
  offset_multiplier_constrain_tests::expect_matvar(x1, mu, sigma);
  offset_multiplier_constrain_tests::expect_matvar(x2, mu, sigma);
}

// real[], real[], real[]
// real[], real, real[]
// real[], real[], real
TEST(mathMixMatFun, offset_multiplier_constrain_matvar_stdvec_constrain) {
  std::vector<double> A{5.0, 2.0, 4.0, -2.0};
  std::vector<double> mum{-3.0, 3.0, -6.0, 6.0};
  std::vector<double> sigmam{-1.0, 5.0, 0.0, 38.0};
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mum, sigmam);
  double mud = -6.0;
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mud, sigmam);
  double sigmad = 8.0;
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mud, sigmad);
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mum, sigmad);
}

// array matrix[], array matrix[], array matrix[]
// array matrix[], array matrix[], matrix[]
// array matrix[], matrix[], array matrix[]
// array matrix[], matrix[], matrix[]
// array matrix[], array matrix[], real
// array matrix[], real, array matrix[]
// array matrix[], matrix[], real
// array matrix[], real, matrix[]
// array matrix[], real, real
TEST(mathMixMatFun, offset_multiplier_matvar_stdvec_mat_scalar_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  // swapping 0.0000001 for 0 causes a failure for the hessian?
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0000001, 0.1;
  Eigen::MatrixXd mu_inner(2, 3);
  mu_inner << -1.0, 1.0, -6.0, 1.0, 0.0, 0.01;
  Eigen::MatrixXd sigma_inner(2, 3);
  sigma_inner << 6.0, 3.0, 12.0, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> mu_vec{mu_inner, mu_inner};
  std::vector<Eigen::MatrixXd> sigma_vec{sigma_inner, sigma_inner};
  double mu_scal = -1.0;
  double sigma_scal = 7.0;
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mu_scal, sigma_inner);
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mu_inner, sigma_scal);
  offset_multiplier_constrain_tests::expect_vec_matvar(A, mu_scal, sigma_scal);
}
