#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(primScalFun, grad_2F2) {
  using stan::math::grad_pFq;

  Eigen::VectorXd a(2);
  a << 4, 2;
  Eigen::VectorXd b(2);
  b << 6, 3;
  double z = 4;

  Eigen::VectorXd grad_a(2);
  Eigen::VectorXd grad_b(2);
  double grad_z;

  grad_pFq(grad_a, grad_b, grad_z, a, b, z);

  EXPECT_FLOAT_EQ(3.924636646666071, grad_a[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, grad_a[1]);

  EXPECT_FLOAT_EQ(-2.775051002566842, grad_b[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, grad_b[1]);

  EXPECT_FLOAT_EQ(4.916522138006060, grad_z);

  grad_a.setZero();

  stan::math::grad_pFq_p(grad_a, a, b, z);

  EXPECT_FLOAT_EQ(3.924636646666071, grad_a[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, grad_a[1]);

  grad_b.setZero();

  stan::math::grad_pFq_q(grad_b, a, b, z);

  EXPECT_FLOAT_EQ(-2.775051002566842, grad_b[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, grad_b[1]);

  grad_z = 0;

  stan::math::grad_pFq_z(grad_z, a, b, z);

  EXPECT_FLOAT_EQ(4.916522138006060, grad_z);

  grad_a.setZero();
  grad_b.setZero();

  stan::math::grad_pFq_pq(grad_a, grad_b, a, b, z);

  EXPECT_FLOAT_EQ(3.924636646666071, grad_a[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, grad_a[1]);
  EXPECT_FLOAT_EQ(-2.775051002566842, grad_b[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, grad_b[1]);

  grad_a.setZero();
  grad_z = 0;
  stan::math::grad_pFq_pz(grad_a, grad_z, a, b, z);
  EXPECT_FLOAT_EQ(3.924636646666071, grad_a[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, grad_a[1]);
  EXPECT_FLOAT_EQ(4.916522138006060, grad_z);

  grad_b.setZero();
  grad_z = 0;
  stan::math::grad_pFq_qz(grad_b, grad_z, a, b, z);
  EXPECT_FLOAT_EQ(-2.775051002566842, grad_b[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, grad_b[1]);
  EXPECT_FLOAT_EQ(4.916522138006060, grad_z);
}

TEST(primScalFun, grad_2F3) {
  using stan::math::grad_pFq;

  Eigen::VectorXd a(2);
  a << 2, 3;
  Eigen::VectorXd b(3);
  b << 2, 4, 5;
  double z = 1;

  Eigen::VectorXd grad_a(2);
  Eigen::VectorXd grad_b(3);
  double grad_z;

  grad_pFq(grad_a, grad_b, grad_z, a, b, z);

  EXPECT_FLOAT_EQ(0.08377717301140296, grad_a[0]);
  EXPECT_FLOAT_EQ(0.05615450733193106, grad_a[1]);

  EXPECT_FLOAT_EQ(-0.08377717301140296, grad_b[0]);
  EXPECT_FLOAT_EQ(-0.04225296806591615, grad_b[1]);
  EXPECT_FLOAT_EQ(-0.03387575989873739, grad_b[2]);

  EXPECT_FLOAT_EQ(0.1712340452215524, grad_z);
}

TEST(primScalFun, grad_4F3) {
  using stan::math::grad_pFq;

  Eigen::VectorXd a(4);
  a << 1, 2, 3, 4;
  Eigen::VectorXd b(3);
  b << 5, 6, 7;
  double z = 1;

  Eigen::VectorXd grad_a(4);
  Eigen::VectorXd grad_b(3);
  double grad_z;

  grad_pFq(grad_a, grad_b, grad_z, a, b, z);

  EXPECT_FLOAT_EQ(0.1587098625610631, grad_a[0]);
  EXPECT_FLOAT_EQ(0.08249338029396848, grad_a[1]);
  EXPECT_FLOAT_EQ(0.05611368752226367, grad_a[2]);
  EXPECT_FLOAT_EQ(0.04261209968272329, grad_a[3]);

  EXPECT_FLOAT_EQ(-0.03438035893346993, grad_b[0]);
  EXPECT_FLOAT_EQ(-0.02882791253333995, grad_b[1]);
  EXPECT_FLOAT_EQ(-0.02482622713079612, grad_b[2]);

  EXPECT_FLOAT_EQ(0.1800529055890911, grad_z);
}

TEST(primScalFun, grad_2F1_derivs_match) {
  using stan::math::grad_pFq_pq;
  using stan::math::grad_2F1;

  Eigen::VectorXd a(2);
  a << 1, 1;
  Eigen::VectorXd b(1);
  b << 1;
  double z = 0.6;

  Eigen::VectorXd grad_a(2);
  Eigen::VectorXd grad_b(1);
  double g_a1;
  double g_b1;

  grad_2F1(g_a1, g_b1, a[0], a[1], b[0], z);
  grad_pFq_pq(grad_a, grad_b, a, b, z);

  EXPECT_FLOAT_EQ(g_a1, grad_a[0]);
  EXPECT_FLOAT_EQ(g_b1, grad_b[0]);
}
