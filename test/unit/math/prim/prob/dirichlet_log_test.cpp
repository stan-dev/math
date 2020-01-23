#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDirichlet, log_matches_lpdf) {
  using stan::math::vector_d;
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha(3, 1);
  alpha << 1.0, 1.0, 1.0;
  EXPECT_FLOAT_EQ(
      (stan::math::dirichlet_lpdf<true, vector_d, vector_d>(theta, alpha)),
      (stan::math::dirichlet_log<true, vector_d, vector_d>(theta, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::dirichlet_lpdf<false, vector_d, vector_d>(theta, alpha)),
      (stan::math::dirichlet_log<false, vector_d, vector_d>(theta, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::dirichlet_lpdf<vector_d, vector_d>(theta, alpha)),
      (stan::math::dirichlet_log<vector_d, vector_d>(theta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpdf(theta, alpha)),
                  (stan::math::dirichlet_log(theta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpdf<true>(theta, alpha)),
                  (stan::math::dirichlet_log<true>(theta, alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpdf<false>(theta, alpha)),
                  (stan::math::dirichlet_log<false>(theta, alpha)));
}
