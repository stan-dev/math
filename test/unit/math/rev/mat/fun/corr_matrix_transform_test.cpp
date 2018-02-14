#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(prob_transform, corr_matrix_jacobian) {
  using stan::math::determinant;
  using stan::math::var;
  using std::fabs;
  using std::log;

  int K = 4;
  int K_choose_2 = 6;
  Matrix<var, Dynamic, 1> X(K_choose_2);
  X << 1.0, 2.0, -3.0, 1.7, 9.8, -1.2;
  std::vector<var> x;
  for (int i = 0; i < X.size(); ++i)
    x.push_back(X(i));
  var lp = 0.0;
  Matrix<var, Dynamic, Dynamic> Sigma
      = stan::math::corr_matrix_constrain(X, K, lp);
  std::vector<var> y;
  for (int m = 0; m < K; ++m)
    for (int n = 0; n < m; ++n)
      y.push_back(Sigma(m, n));
  EXPECT_EQ(K_choose_2, y.size());

  std::vector<std::vector<double> > j;
  stan::math::jacobian(y, x, j);

  Matrix<double, Dynamic, Dynamic> J(X.size(), X.size());
  for (int m = 0; m < J.rows(); ++m)
    for (int n = 0; n < J.cols(); ++n)
      J(m, n) = j[m][n];

  double log_abs_jacobian_det = log(fabs(determinant(J)));
  EXPECT_FLOAT_EQ(log_abs_jacobian_det, lp.val());
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  int K = 4;
  int K_choose_2 = 6;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> X(K_choose_2);
  X << 1.0, 2.0, -3.0, 1.7, 9.8, -1.2;
  stan::math::var lp = 0.0;
  test::check_varis_on_stack(stan::math::corr_matrix_constrain(X, K, lp));
  test::check_varis_on_stack(stan::math::corr_matrix_constrain(X, K));
}
