#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;
using stan::math::multi_student_t_log;

TEST(ProbDistributionsMultiStudentT,fvar_var) {
  using stan::math::var;
  using stan::math::fvar;
  Matrix<fvar<var>,Dynamic,1> y(3,1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<var>,Dynamic,1> mu(3,1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<var>,Dynamic,Dynamic> Sigma(3,3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  double nu = 4.0;

  for (int i = 0; i < 3; i++) {
    y(i).d_ = 1.0;
    mu(i).d_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i,j).d_ = 1.0;
  }

  fvar<var> lp = multi_student_t_log(y,nu,mu,Sigma);
  EXPECT_NEAR(-10.1246,lp.val_.val(),0.0001);
  EXPECT_NEAR(-0.0411685,lp.d_.val(),0.0001);
}

TEST(ProbDistributionsMultiStudentT,fvar_fvar_var) {
  using stan::math::var;
  using stan::math::fvar;
  Matrix<fvar<fvar<var> >,Dynamic,1> y(3,1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<fvar<var> >,Dynamic,1> mu(3,1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<fvar<var> >,Dynamic,Dynamic> Sigma(3,3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  double nu = 4.0;

  for (int i = 0; i < 3; i++) {
    y(i).d_.val_ = 1.0;
    mu(i).d_.val_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i,j).d_.val_ = 1.0;
  }

  fvar<fvar<var> > lp = multi_student_t_log(y,nu,mu,Sigma);
  EXPECT_NEAR(-10.1246,lp.val_.val_.val(),0.0001);
  EXPECT_NEAR(-0.0411685,lp.d_.val_.val(),0.0001);
}
