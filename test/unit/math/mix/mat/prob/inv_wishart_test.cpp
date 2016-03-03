#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>



#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/digamma.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;

using stan::math::inv_wishart_log;

TEST(ProbDistributionsInvWishart,fvar_var) {
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<var>,Dynamic,Dynamic> Y(3,3);
  Y <<  12.147233, -11.9036079, 1.0910458,
    -11.903608,  16.7585782, 0.8530256,
    1.091046,   0.8530256, 2.5786609;

  Matrix<fvar<var>,Dynamic,Dynamic> Sigma(3,3);
  Sigma << 7.785215,  3.0597878,  1.1071663,
    3.059788, 10.3515035, -0.1232598,
    1.107166, -0.1232598,  7.7623386;
  
  double dof = 4.0;
  double log_p = log(2.008407e-08);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Y(i,j).d_ = 1.0;
      Sigma(i,j).d_ = 1.0;
    }

  EXPECT_NEAR(log_p, stan::math::inv_wishart_log(Y,dof,Sigma).val_.val(), 0.01);
  EXPECT_NEAR(-1.4893348387330674, stan::math::inv_wishart_log(Y,dof,Sigma).d_.val(), 0.01);
}

TEST(ProbDistributionsInvWishart,fvar_fvar_var) {
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<fvar<var> >,Dynamic,Dynamic> Y(3,3);
  Y <<  12.147233, -11.9036079, 1.0910458,
    -11.903608,  16.7585782, 0.8530256,
    1.091046,   0.8530256, 2.5786609;

  Matrix<fvar<fvar<var> >,Dynamic,Dynamic> Sigma(3,3);
  Sigma << 7.785215,  3.0597878,  1.1071663,
    3.059788, 10.3515035, -0.1232598,
    1.107166, -0.1232598,  7.7623386;
  
  double dof = 4.0;
  double log_p = log(2.008407e-08);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Y(i,j).d_ = 1.0;
      Sigma(i,j).d_ = 1.0;
    }

  EXPECT_NEAR(log_p, stan::math::inv_wishart_log(Y,dof,Sigma).val_.val_.val(), 0.01);
  EXPECT_NEAR(-1.4893348387330674, stan::math::inv_wishart_log(Y,dof,Sigma).d_.val_.val(), 0.01);
}
