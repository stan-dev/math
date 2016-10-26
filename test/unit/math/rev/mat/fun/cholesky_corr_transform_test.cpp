#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>

void 
test_cholesky_correlation_jacobian(const Eigen::Matrix<stan::math::var,
                                                       Eigen::Dynamic,1>& y,
                                   int K) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;
  using stan::math::cholesky_corr_constrain;

  int K_choose_2 = (K * (K - 1)) / 2;

  vector<var> indeps;
  for (int i = 0; i < y.size(); ++i)
    indeps.push_back(y(i));

  var lp = 0;
  Matrix<var,Dynamic,Dynamic> x
    = cholesky_corr_constrain(y,K,lp);

  vector<var> deps;
  for (int i = 1; i < K; ++i)
    for (int j = 0; j < i; ++j)
      deps.push_back(x(i,j));
  
  vector<vector<double> > jacobian;
  stan::math::jacobian(deps,indeps,jacobian);

  Matrix<double,Dynamic,Dynamic> J(K_choose_2,K_choose_2);
  for (int m = 0; m < K_choose_2; ++m)
    for (int n = 0; n < K_choose_2; ++n)
      J(m,n) = jacobian[m][n];

  
  double det_J = J.determinant();
  double log_det_J = log(fabs(det_J));

  EXPECT_FLOAT_EQ(log_det_J, lp.val()) << "J = " << J << std::endl << "det_J = " << det_J;
  
}

TEST(probTransform,choleskyCorrJacobian) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;

  // K = 1; (K choose 2) = 0
  Matrix<var,Dynamic,1> y1;
  EXPECT_EQ(0,y1.size());
  test_cholesky_correlation_jacobian(y1,1);

  // K = 2; (K choose 2) = 1
  Matrix<var,Dynamic,1> y2(1);
  y2 << -1.7;
  test_cholesky_correlation_jacobian(y2,2);

  // K = 3; (K choose 2) = 3
  Matrix<var,Dynamic,1> y3(3);
  y3 << -1.7, 2.9, 0.01;
  test_cholesky_correlation_jacobian(y3,3);

  // K = 4;  (K choose 2) = 6
  Matrix<var,Dynamic,1> y4(6);
  y4 << 1.0, 2.0, -3.0, 1.5, 0.2, 2.0;
  test_cholesky_correlation_jacobian(y4,4);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::vector_v y(3);
  y << -1.7, 2.9, 0.01;
  stan::math::var lp(0);

  test::check_varis_on_stack(stan::math::cholesky_corr_constrain(y, 3, lp));
  test::check_varis_on_stack(stan::math::cholesky_corr_constrain(y, 3));
}
