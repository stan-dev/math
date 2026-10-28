#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, inc_beta_dda_vectorized) {
  using stan::math::digamma;
  using stan::math::inc_beta_dda;

  Eigen::VectorXd a(4);
  a << 1.145, 446, 2, 600;
  Eigen::VectorXd b(4);
  b << 6786, 5, 3, 600;
  Eigen::VectorXd z(4);
  z << 1.8e-4, 0.76, 0.25, 0.5;
  Eigen::VectorXd digamma_a(4);
  Eigen::VectorXd digamma_ab(4);
  for (int i = 0; i < 4; ++i) {
    digamma_a[i] = digamma(a[i]);
    digamma_ab[i] = digamma(a[i] + b[i]);
  }

  Eigen::ArrayXd all_eigen = inc_beta_dda(a, b, z, digamma_a, digamma_ab);
  ASSERT_EQ(4, all_eigen.size());
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(inc_beta_dda(a[i], b[i], z[i], digamma_a[i], digamma_ab[i]),
              all_eigen[i]);
  }

  std::vector<double> z_std{1.8e-4, 0.76, 0.25, 0.5};
  Eigen::ArrayXd std_z = inc_beta_dda(a, b, z_std, digamma_a, digamma_ab);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(all_eigen[i], std_z[i]);
  }

  // scalar shapes broadcast against a vector of upper limits
  Eigen::ArrayXd only_z = inc_beta_dda(2.0, 3.0, z, digamma(2.0), digamma(5.0));
  ASSERT_EQ(4, only_z.size());
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(inc_beta_dda(2.0, 3.0, z[i], digamma(2.0), digamma(5.0)),
              only_z[i]);
  }

  Eigen::VectorXd too_short(3);
  too_short << 0.1, 0.2, 0.3;
  EXPECT_THROW(inc_beta_dda(a, b, too_short, digamma_a, digamma_ab),
               std::invalid_argument);
}
