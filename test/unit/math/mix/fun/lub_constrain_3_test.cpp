#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/lub_constrain_helpers.hpp>

TEST(mathMixMatFun, lub_constrain_vector_scalar_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  double lb = -2.0;
  Eigen::MatrixXd ub(2, 2);
  ub << -1.0, 5.0, 0.0, 38.0;

  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // ub inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // lb inf
  ub(1, 1) = 38.0;
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);
}

TEST(mathMixMatFun, lub_constrain_vector_vector_vector) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.0000001, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 0.0, -6.0, 6.0;
  Eigen::MatrixXd ub(2, 2);
  ub << -1.0, 5.0, 0.0, 38.0;
  // lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // ub inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // lb inf
  ub(1, 1) = 38.0;
  lb(1, 1) = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // both inf
  ub(1, 1) = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
}
