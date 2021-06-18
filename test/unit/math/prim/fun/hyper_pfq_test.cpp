#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(primScalFun, hyper_pfq) {
  using stan::math::grad_pFq;

  Eigen::VectorXd ab(2);
  ab << 4,2;
  Eigen::VectorXd cd(2);
  cd << 6,3;
  double z = 4;
  auto t = grad_pFq(ab, cd, z);

  Eigen::VectorXd d_p = std::get<0>(t);
  Eigen::VectorXd d_q = std::get<1>(t);
  double d_z = std::get<2>(t);

  EXPECT_FLOAT_EQ(3.924636646666071, d_p[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, d_p[1]);
  
  EXPECT_FLOAT_EQ(-2.775051002566842, d_q[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, d_q[1]);

  EXPECT_FLOAT_EQ(4.916522138006060, d_z);
}
