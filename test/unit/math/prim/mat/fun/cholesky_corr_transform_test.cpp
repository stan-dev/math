#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(ProbTransform,CholeskyCorrelation4) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  Matrix<double,Dynamic,Dynamic> L(4,4);
  L << 
    1, 0, 0, 0,
    -0.2, 0.9797959, 0, 0,
    0.5, -0.3, 0.8124038, 0,
    0.7, -0.2, 0.6, 0.3316625;

  Matrix<double,Dynamic,1> y
    = stan::math::cholesky_corr_free(L);

  Matrix<double,Dynamic,Dynamic> x
    = stan::math::cholesky_corr_constrain(y,4);
  
  Matrix<double,Dynamic,1> yrt
    = stan::math::cholesky_corr_free(x);

  EXPECT_EQ(y.size(), yrt.size());
  for (int i = 0; i < yrt.size(); ++i)
    EXPECT_FLOAT_EQ(y(i), yrt(i));

  for (int m = 0; m < 4; ++m)
    for (int n = 0; n < 4; ++n)
      EXPECT_FLOAT_EQ(L(m,n), x(m,n));
}

void 
test_cholesky_correlation_values(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& L) {
  using std::vector;
  using stan::math::cholesky_corr_constrain;
  using stan::math::cholesky_corr_free;
  int K = L.rows();
  int K_choose_2 = (K * (K - 1)) / 2;

  // test number of free parameters
  Matrix<double,Dynamic,1> y
    = stan::math::cholesky_corr_free(L);
  EXPECT_EQ(K_choose_2, y.size());

  // test transform roundtrip without Jacobian
  Matrix<double,Dynamic,Dynamic> x
    = stan::math::cholesky_corr_constrain(y,K);
  
  Matrix<double,Dynamic,1> yrt
    = stan::math::cholesky_corr_free(x);

  EXPECT_EQ(y.size(), yrt.size());
  for (int i = 0; i < yrt.size(); ++i)
    EXPECT_FLOAT_EQ(y(i), yrt(i));

  for (int m = 0; m < K; ++m)
    for (int n = 0; n < K; ++n)
      EXPECT_FLOAT_EQ(L(m,n), x(m,n));


  // test transform roundtrip with Jacobian (Jacobian itself tested above)
  double lp;
  Matrix<double,Dynamic,Dynamic> x2
    = stan::math::cholesky_corr_constrain(y,K,lp);
  
  Matrix<double,Dynamic,1> yrt2
    = stan::math::cholesky_corr_free(x2);

  EXPECT_EQ(y.size(), yrt2.size());
  for (int i = 0; i < yrt2.size(); ++i)
    EXPECT_FLOAT_EQ(y(i), yrt2(i));

  for (int m = 0; m < K; ++m)
    for (int n = 0; n < K; ++n)
      EXPECT_FLOAT_EQ(L(m,n), x2(m,n));
}

TEST(ProbTransform,CholeskyCorrelationRoundTrips) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<double,Dynamic,Dynamic> L1(1,1);
  L1 << 1;
  test_cholesky_correlation_values(L1);

  Matrix<double,Dynamic,Dynamic> L2(2,2);
  L2 << 
    1, 0,
    -0.5, 0.8660254;
  test_cholesky_correlation_values(L2);
    
  Matrix<double,Dynamic,Dynamic> L4(4,4);
  L4 << 
    1, 0, 0, 0,
    -0.2, 0.9797959, 0, 0,
    0.5, -0.3, 0.8124038, 0,
    0.7, -0.2, 0.6, 0.3316625;
  test_cholesky_correlation_values(L4);
}
