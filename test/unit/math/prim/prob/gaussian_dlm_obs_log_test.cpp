#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbGaussianDlmObs, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y(1, 10);
  y << -0.286804393606091, 1.30654039013044, 0.184631538931975,
      1.76116251447979, 1.64691178557684, 0.0599998209370169,
      -0.498099220647035, 1.77794756092381, -0.435458550812876,
      1.17332931763075;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F(1, 1);
  F << 0.585528817843856;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> G(1, 1);
  G << -0.109303314681054;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> V(1, 1);
  V << 2.25500747900521;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> W(1, 1);
  W << 0.461487989960454;
  Eigen::Matrix<double, Eigen::Dynamic, 1> m0(1);
  m0 << 11.5829455171551;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C0(1, 1);
  C0 << 65.2373490156606;

  EXPECT_FLOAT_EQ(
      (stan::math::gaussian_dlm_obs_lpdf<true>(y, F, G, V, W, m0, C0)),
      (stan::math::gaussian_dlm_obs_log<true>(y, F, G, V, W, m0, C0)));
  EXPECT_FLOAT_EQ(
      (stan::math::gaussian_dlm_obs_lpdf<false>(y, F, G, V, W, m0, C0)),
      (stan::math::gaussian_dlm_obs_log<false>(y, F, G, V, W, m0, C0)));
  EXPECT_FLOAT_EQ((stan::math::gaussian_dlm_obs_lpdf<>(y, F, G, V, W, m0, C0)),
                  (stan::math::gaussian_dlm_obs_log<>(y, F, G, V, W, m0, C0)));
  EXPECT_FLOAT_EQ((stan::math::gaussian_dlm_obs_lpdf(y, F, G, V, W, m0, C0)),
                  (stan::math::gaussian_dlm_obs_log(y, F, G, V, W, m0, C0)));
  EXPECT_FLOAT_EQ(
      (stan::math::gaussian_dlm_obs_lpdf<true>(y, F, G, V, W, m0, C0)),
      (stan::math::gaussian_dlm_obs_log<true>(y, F, G, V, W, m0, C0)));
  EXPECT_FLOAT_EQ(
      (stan::math::gaussian_dlm_obs_lpdf<false>(y, F, G, V, W, m0, C0)),
      (stan::math::gaussian_dlm_obs_log<false>(y, F, G, V, W, m0, C0)));
}
