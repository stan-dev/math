#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(prob_transform,ordered_jacobian_ad) {
  using stan::math::var;
  using stan::math::ordered_constrain;
  using stan::math::determinant;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<double,Dynamic,1> x(3);
  x << -12.0, 3.0, -1.9;
  double lp = 0.0;
  Matrix<double,Dynamic,1> y = ordered_constrain(x,lp);

  Matrix<var,Dynamic,1> xv(3);
  xv << -12.0, 3.0, -1.9;

  std::vector<var> xvec(3);
  for (int i = 0; i < 3; ++i)
    xvec[i] = xv[i];

  Matrix<var,Dynamic,1> yv = ordered_constrain(xv);


  EXPECT_EQ(y.size(), yv.size());
  for (int i = 0; i < y.size(); ++i)
    EXPECT_FLOAT_EQ(y(i),yv(i).val());

  std::vector<var> yvec(3);
  for (unsigned int i = 0; i < 3; ++i)
    yvec[i] = yv[i];

  std::vector<std::vector<double> > j;
  stan::math::jacobian(yvec,xvec,j);

  Matrix<double,Dynamic,Dynamic> J(3,3);
  for (int m = 0; m < 3; ++m)
    for (int n = 0; n < 3; ++n)
      J(m,n) = j[m][n];
  
  double log_abs_jacobian_det = log(fabs(determinant(J)));
  EXPECT_FLOAT_EQ(log_abs_jacobian_det, lp);
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  Eigen::Matrix<stan::math::var,Eigen::Dynamic,1> x(3);

  x << -12.0, 3.0, -1.9;
  stan::math::var lp = 0.0;

  test::check_varis_on_stack(stan::math::ordered_constrain(x, lp));
  test::check_varis_on_stack(stan::math::ordered_constrain(x));
}
