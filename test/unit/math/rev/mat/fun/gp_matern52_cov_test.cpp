#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <stan/math/rev/mat.hpp>
#include <string>
#include <vector>
#include <type_traits>

TEST(MathPrimMat, vec_double_gp_matern52_cov1) {
  stan::math::var sigma = 0.2;
  stan::math::var l = 5.0;

  std::vector<stan::math::var> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          stan::math::value_of(
              sigma * sigma
              * (1
                 + std::pow(5, 0.5) / l
                       * stan::math::sqrt(
                             stan::math::squared_distance(x[i], x[j]))
                 + (5.0 / 3.0) * stan::math::squared_distance(x[i], x[j])
                       / std::pow(value_of(l), 2))
              * std::exp(stan::math::value_of(
                    -1.0 * pow(5.0, 0.5)
                    * stan::math::sqrt(stan::math::squared_distance(x[i], x[j]))
                    / l))),
          cov(i, j).val())
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, output_type_checking) {
  stan::math::var sigma = 0.2;
  stan::math::var l = 5.0;

  std::vector<stan::math::var> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov1;
  EXPECT_NO_THROW(cov1 = stan::math::gp_matern52_cov(value_of(x),
                                                     value_of(sigma),
                                                     value_of(l)));
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov2;
  EXPECT_NO_THROW(cov2 = stan::math::gp_matern52_cov(x, value_of(sigma),
                                                     value_of(l)));
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov3;
  EXPECT_NO_THROW(cov3 = stan::math::gp_matern52_cov(value_of(x), sigma,
                                                     value_of(l)));
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov4;
  EXPECT_NO_THROW(cov3 = stan::math::gp_matern52_cov(value_of(x),
                                                     value_of(sigma), l));
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_TRUE((std::is_same<double, decltype(cov1(i, j).val())>::value));
      // EXPECT_TRUE((boost::is_same<stan::math::var, decltype(cov2(i, j))>::value));
      //      EXPECT_TRUE((std::is_same<stan::math::var, decltype(value_of(cov2(i, j)))>::value));
    }
  }
}

