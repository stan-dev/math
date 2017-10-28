#include <stan/math/prim/mat/fun/LDLT_factor.hpp>
#include <stan/math/prim/mat/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrix, LDLTfactor_Values) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> spd(3, 3);
  spd << 2.0, -1.0, 0.0,
        -1.0, 2.0, -1.0,
         0.0, -1.0, 2.0;

  stan::math::LDLT_factor<double,
                          Eigen::Dynamic, Eigen::Dynamic> stan_ldlt_1(spd);
  stan::math::LDLT_factor<double, Eigen::Dynamic, Eigen::Dynamic> stan_ldlt_2;
  stan_ldlt_2.compute(spd);

  Eigen::LDLT<Eigen::Matrix<double,
                            Eigen::Dynamic, Eigen::Dynamic>> eigen_ldlt(spd);

  double stan_ldlt_1_logdet = stan_ldlt_1.log_abs_det();
  double stan_ldlt_2_logdet = stan_ldlt_2.log_abs_det();
  double eigen_ldlt_logdet = eigen_ldlt.vectorD().array().log().sum();

  EXPECT_EQ(stan_ldlt_1_logdet, eigen_ldlt_logdet);
  EXPECT_EQ(stan_ldlt_1_logdet, stan_ldlt_2_logdet);
}

TEST(MathMatrix, LDLTfactor_InvalidValues) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double infty = std::numeric_limits<double>::infinity();
  double neg_infty = -std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> non_spd(3, 3);
  non_spd << 2.0, 0.0, 0.0,
            -1.0, 0.0, -1.0,
             0.0, -1.0, 2.0;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> spd_nan(3, 3);
  spd_nan << 2.0, -1.0, nan,
            -1.0, 2.0, -1.0,
             nan, -1.0, 2.0;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> spd_infty(3, 3);
  spd_infty << 2.0, -1.0, infty,
              -1.0, 2.0, -1.0,
               infty, -1.0, 2.0;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> spd_neg_infty(3, 3);
  spd_neg_infty << 2.0, -1.0, neg_infty,
                  -1.0, 2.0, -1.0,
                   neg_infty, -1.0, 2.0;

  stan::math::LDLT_factor<double,
                          Eigen::Dynamic,
                          Eigen::Dynamic> non_spd_ldlt(non_spd);
  stan::math::LDLT_factor<double,
                          Eigen::Dynamic,
                          Eigen::Dynamic> nan_ldlt(spd_nan);
  stan::math::LDLT_factor<double,
                          Eigen::Dynamic,
                          Eigen::Dynamic> infty_ldlt(spd_infty);
  stan::math::LDLT_factor<double,
                          Eigen::Dynamic,
                          Eigen::Dynamic> neg_infty_ldlt(spd_neg_infty);

  double non_spd_ldlt_logdet = non_spd_ldlt.log_abs_det();
  double nan_ldlt_logdet = nan_ldlt.log_abs_det();
  double infty_ldlt_logdet = infty_ldlt.log_abs_det();
  double neg_infty_ldlt_logdet = neg_infty_ldlt.log_abs_det();

  EXPECT_TRUE(isnan(non_spd_ldlt_logdet));
  EXPECT_TRUE(isnan(nan_ldlt_logdet));
  EXPECT_TRUE(isnan(infty_ldlt_logdet));
  EXPECT_TRUE(isnan(neg_infty_ldlt_logdet));
}
