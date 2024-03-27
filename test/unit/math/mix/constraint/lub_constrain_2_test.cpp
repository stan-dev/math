#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/constraint/lub_constrain_helpers.hpp>

TEST(mathMixMatFun, lub_constrain_vector_scalar_scalar) {
  Eigen::MatrixXd x1(1, 1);
  x1 << 0.7;
  Eigen::MatrixXd x2(1, 1);
  x2 << -1.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);

  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect(x2, lb_inf, ub_inf);
}

TEST(mathMixMatFun, lub_constrain_vector_vector_scalar) {
  Eigen::MatrixXd x1(2, 2);
  x1 << 5.0, 2.0, 4.0, 5.0;
  Eigen::MatrixXd x2(2, 2);
  x2 << -1.1, 0.005, 1.0, 3.0;
  Eigen::MatrixXd lb(2, 2);
  lb << -3.0, 5.0, -6.0, 6.0;
  double ub = 13.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  lb(1, 1) = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);
}
