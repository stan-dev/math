#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbOrderedLogistic, log_matches_lpmf) {
  /**
   * Variables for int ~ ordered_logistic(double, vector)
   */
  int K = 5;
  Eigen::Matrix<double, Eigen::Dynamic, 1> c(K - 1);
  c << -1.7, -0.3, 1.2, 2.6;
  double lambda = 1.1;

  /**
   * Variables for int[] ~ ordered_logistic(vector, vector)
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> lamvec(K);
  lamvec << 1.1, 1.3, 2.2, 3.1, 2.1;

  std::vector<int> kvec{1, 2, 3, 1, 3};

  /**
   * Variables for int[] ~ ordered_logistic(vector, vector[])
   */

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_1(K - 1);
  c_1 << -2.3, -0.3, 1.5, 2.6;

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_2(K - 1);
  c_2 << -1.7, -0.3, 1.2, 2.6;

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_3(K - 1);
  c_3 << -1.7, 0.2, 1.4, 2.6;

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_4(K - 1);
  c_4 << -1.5, -0.3, 0.2, 2.6;

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_5(K - 1);
  c_5 << -1.7, -0.3, 1.7, 2.3;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> cutvec{c_1, c_2, c_3,
                                                               c_4, c_5};

  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf(3, lambda, c)),
                  (stan::math::ordered_logistic_log(3, lambda, c)));

  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf(kvec, lamvec, c)),
                  (stan::math::ordered_logistic_log(kvec, lamvec, c)));

  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf(kvec, lamvec, cutvec)),
                  (stan::math::ordered_logistic_log(kvec, lamvec, cutvec)));

  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf<true>(3, lambda, c)),
                  (stan::math::ordered_logistic_log<true>(3, lambda, c)));
  EXPECT_FLOAT_EQ(
      (stan::math::ordered_logistic_lpmf<double, double>(3, lambda, c)),
      (stan::math::ordered_logistic_log<double, double>(3, lambda, c)));
}
