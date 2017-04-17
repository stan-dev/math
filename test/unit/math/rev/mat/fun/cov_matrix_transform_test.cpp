#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(prob_transform,cov_matrix_jacobian) {
  using stan::math::var;
  using stan::math::determinant;
  using std::log;
  using std::fabs;

  int K = 4;
  //unsigned int K = 4;
  unsigned int K_choose_2 = 6;
  Matrix<var,Dynamic,1> X(K_choose_2 + K);
  X << 1.0, 2.0, -3.0, 1.7, 9.8, 
    -12.2, 0.4, 0.2, 1.2, 2.7;
  std::vector<var> x;
  for (int i = 0; i < X.size(); ++i)
    x.push_back(X(i));
  var lp = 0.0;
  Matrix<var,Dynamic,Dynamic> Sigma = stan::math::cov_matrix_constrain(X,K,lp);
  std::vector<var> y;
  for (int m = 0; m < K; ++m)
    for (int n = 0; n <= m; ++n)
      y.push_back(Sigma(m,n));

  std::vector<std::vector<double> > j;
  stan::math::jacobian(y,x,j);

  Matrix<double,Dynamic,Dynamic> J(10,10);
  for (int m = 0; m < 10; ++m)
    for (int n = 0; n < 10; ++n)
      J(m,n) = j[m][n];

  double log_abs_jacobian_det = log(fabs(determinant(J)));
  EXPECT_FLOAT_EQ(log_abs_jacobian_det,lp.val());
}

TEST(prob_transform, check_varis_on_stack) {
  using stan::math::var;
  int K = 4;
  unsigned int K_choose_2 = 6;
  Eigen::Matrix<var,Eigen::Dynamic,1> X(K_choose_2 + K);
  X << 1.0, 2.0, -3.0, 1.7, 9.8, 
    -12.2, 0.4, 0.2, 1.2, 2.7;
  var lp = 0.0;

  test::check_varis_on_stack(stan::math::cov_matrix_constrain(X, K, lp));
  test::check_varis_on_stack(stan::math::cov_matrix_constrain(X, K));
}
