#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, inc_beta_ddz) {
  using stan::math::inc_beta_ddz;

  double small_a = 1.5;
  double large_a = 15000;

  double small_b = 1.25;
  double large_b = 12500;

  double small_z = 0.001;
  double mid_z = 0.5;
  double large_z = 0.999;

  EXPECT_FLOAT_EQ(0.063300692, inc_beta_ddz(small_a, small_b, small_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(1.1905416, inc_beta_ddz(small_a, small_b, mid_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.35587692, inc_beta_ddz(small_a, small_b, large_z))
      << "reasonable values for a, b, x";

  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(large_a, small_b, small_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(large_a, small_b, mid_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.009898182, inc_beta_ddz(large_a, small_b, large_z))
      << "reasonable values for a, b, x";

  EXPECT_FLOAT_EQ(0.1848717, inc_beta_ddz(small_a, large_b, small_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(small_a, large_b, mid_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(small_a, large_b, large_z))
      << "reasonable values for a, b, x";

  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(large_a, large_b, small_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(large_a, large_b, mid_z))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddz(large_a, large_b, large_z))
      << "reasonable values for a, b, x";
}

TEST(MathFunctions, inc_beta_ddz_infinite_density) {
  using stan::math::inc_beta_ddz;
  // the density is infinite at z = 1 with b < 1 and at z = 0 with a < 1
  double at_one = 0;
  double at_zero = 0;
  EXPECT_NO_THROW(at_one = inc_beta_ddz(1.05, 0.35, 1.0));
  EXPECT_NO_THROW(at_zero = inc_beta_ddz(0.5, 2.0, 0.0));
  EXPECT_TRUE(std::isinf(at_one) && at_one > 0);
  EXPECT_TRUE(std::isinf(at_zero) && at_zero > 0);
}

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
