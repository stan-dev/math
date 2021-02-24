#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/lub_constrain_helpers.hpp>

// array matrix[], array matrix[], array matrix[]
// array matrix[], array matrix[], matrix[]
// array matrix[], matrix[], array matrix[]
// array matrix[], matrix[], matrix[]
// array matrix[], array matrix[], real
// array matrix[], real, array matrix[]
// array matrix[], matrix[], real
// array matrix[], real, matrix[]
// array matrix[], real, real
TEST(mathMixMatFun, lub_stdvec_mat_scalar1_constrain) {
  Eigen::MatrixXd A_inner(2, 3);
  // swapping 0.0000001 for 0 causes a failure for the hessian?
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.0000001, 0.1;
  Eigen::MatrixXd lb_inner(2, 3);
  lb_inner << -1.0, 1.0, -6.0, 1.0, 0.0, 0.01;
  Eigen::MatrixXd ub_inner(2, 3);
  ub_inner << 6.0, 3.0, 12.0, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> lb_vec{lb_inner, lb_inner};
  std::vector<Eigen::MatrixXd> ub_vec{ub_inner, ub_inner};
  double lb_scal = -1.0;
  double ub_scal = 7.0;
  lub_constrain_tests::expect_vec(A, lb_vec, ub_scal);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_vec);
}

TEST(mathMixMatFun, lub_stdvec_mat_scalar1_constrain_infty) {
  Eigen::MatrixXd A_inner(2, 3);
  A_inner << 5.0, 2.0, 4.0, -2.0, 0.05, 0.1;
  Eigen::MatrixXd lb_inner(2, 3);
  lb_inner << stan::math::NEGATIVE_INFTY, 1.0, stan::math::NEGATIVE_INFTY, 1.0,
      0.0, 0.01;
  Eigen::MatrixXd ub_inner(2, 3);
  ub_inner << 6.0, stan::math::INFTY, stan::math::INFTY, 38.0, 0.1, 0.15;

  std::vector<Eigen::MatrixXd> A{A_inner, A_inner};
  std::vector<Eigen::MatrixXd> lb_vec{lb_inner, lb_inner};
  std::vector<Eigen::MatrixXd> ub_vec{ub_inner, ub_inner};
  double lb_scal = stan::math::NEGATIVE_INFTY;
  double ub_scal = stan::math::INFTY;
  lub_constrain_tests::expect_vec(A, lb_vec, ub_scal);
  lub_constrain_tests::expect_vec(A, lb_scal, ub_vec);
}
