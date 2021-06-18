#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(primScalFun, grad_2F2) {
  using stan::math::grad_pFq;

  Eigen::VectorXd p(2);
  p << 4, 2;
  Eigen::VectorXd q(2);
  q << 6, 3;
  double z = 4;

  Eigen::VectorXd grad_p(2);
  Eigen::VectorXd grad_q(2);
  double grad_z;

  grad_pFq(grad_p, grad_q, grad_z, p, q, z);

  EXPECT_FLOAT_EQ(3.924636646666071, grad_p[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, grad_p[1]);
  
  EXPECT_FLOAT_EQ(-2.775051002566842, grad_q[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, grad_q[1]);

  EXPECT_FLOAT_EQ(4.916522138006060, grad_z);
}

TEST(primScalFun, grad_2F3) {
  using stan::math::grad_pFq;

  Eigen::VectorXd p(2);
  p << 2, 3;
  Eigen::VectorXd q(3);
  q << 2, 4, 5;
  double z = 1;
  
  Eigen::VectorXd grad_p(2);
  Eigen::VectorXd grad_q(3);
  double grad_z;
  
  grad_pFq(grad_p, grad_q, grad_z, p, q, z);

  EXPECT_FLOAT_EQ(0.08377717301140296, grad_p[0]);
  EXPECT_FLOAT_EQ(0.05615450733193106, grad_p[1]);
  
  EXPECT_FLOAT_EQ(-0.08377717301140296, grad_q[0]);
  EXPECT_FLOAT_EQ(-0.04225296806591615, grad_q[1]);
  EXPECT_FLOAT_EQ(-0.03387575989873739, grad_q[2]);

  EXPECT_FLOAT_EQ(0.1712340452215524, grad_z);
}

TEST(primScalFun, grad_4F3) {
  using stan::math::grad_pFq;

  Eigen::VectorXd p(4);
  p << 1, 2, 3, 4;
  Eigen::VectorXd q(3);
  q << 5, 6, 7;
  double z = 1;

  Eigen::VectorXd grad_p(4);
  Eigen::VectorXd grad_q(3);
  double grad_z;

  grad_pFq(grad_p, grad_q, grad_z, p, q, z);

  EXPECT_FLOAT_EQ(0.1587098625610631, grad_p[0]);
  EXPECT_FLOAT_EQ(0.08249338029396848, grad_p[1]);
  EXPECT_FLOAT_EQ(0.05611368752226367, grad_p[2]);
  EXPECT_FLOAT_EQ(0.04261209968272329, grad_p[3]);
  
  EXPECT_FLOAT_EQ(-0.03438035893346993, grad_q[0]);
  EXPECT_FLOAT_EQ(-0.02882791253333995, grad_q[1]);
  EXPECT_FLOAT_EQ(-0.02482622713079612, grad_q[2]);

  EXPECT_FLOAT_EQ(0.1800529055890911, grad_z);
}