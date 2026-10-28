#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, inc_beta_ddz_vectorized) {
  using stan::math::inc_beta_ddz;

  Eigen::VectorXd a(3);
  a << 1.5, 15.0, 0.75;
  Eigen::VectorXd b(3);
  b << 1.25, 4.0, 30.0;
  Eigen::VectorXd z(3);
  z << 0.001, 0.5, 0.999;

  Eigen::VectorXd all_eigen = inc_beta_ddz(a, b, z);
  ASSERT_EQ(3, all_eigen.size());
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(inc_beta_ddz(a[i], b[i], z[i]), all_eigen[i]);
  }

  std::vector<double> a_std{1.5, 15.0, 0.75};
  std::vector<double> b_std{1.25, 4.0, 30.0};
  std::vector<double> z_std{0.001, 0.5, 0.999};
  std::vector<double> all_std = inc_beta_ddz(a_std, b_std, z_std);
  ASSERT_EQ(3, all_std.size());
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(inc_beta_ddz(a_std[i], b_std[i], z_std[i]), all_std[i]);
  }

  Eigen::VectorXd scalar_a = inc_beta_ddz(2.5, b, z);
  Eigen::VectorXd scalar_b = inc_beta_ddz(a, 3.5, z);
  Eigen::VectorXd scalar_z = inc_beta_ddz(a, b, 0.3);
  Eigen::VectorXd only_z = inc_beta_ddz(2.5, 3.5, z);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(inc_beta_ddz(2.5, b[i], z[i]), scalar_a[i]);
    EXPECT_EQ(inc_beta_ddz(a[i], 3.5, z[i]), scalar_b[i]);
    EXPECT_EQ(inc_beta_ddz(a[i], b[i], 0.3), scalar_z[i]);
    EXPECT_EQ(inc_beta_ddz(2.5, 3.5, z[i]), only_z[i]);
  }
}
