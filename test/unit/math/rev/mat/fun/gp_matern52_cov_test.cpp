#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>
#include <cmath>
#include <limits>
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

  // If the scalar types of all arguments are doubles, the scalar type of the
  // output should be a double
  EXPECT_TRUE(
      (std::is_same<Eigen::MatrixXd,
                    decltype(stan::math::gp_matern52_cov(
                        value_of(x), value_of(sigma), value_of(l)))>::value));

  // If the scalar types of any argument is a var, the scalar type of the output
  // should be a var
  EXPECT_TRUE((std::is_same<
               Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>,
               decltype(stan::math::gp_matern52_cov(x, sigma,
                                                    value_of(l)))>::value));
  EXPECT_TRUE((std::is_same<
               Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>,
               decltype(stan::math::gp_matern52_cov(x, value_of(sigma),
                                                    l))>::value));
  EXPECT_TRUE((std::is_same<
               Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>,
               decltype(stan::math::gp_matern52_cov(x, value_of(sigma),
                                                    value_of(l)))>::value));
  EXPECT_TRUE((std::is_same<
               Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>,
               decltype(stan::math::gp_matern52_cov(value_of(x), sigma,
                                                    l))>::value));
  EXPECT_TRUE((std::is_same<
               Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>,
               decltype(stan::math::gp_matern52_cov(value_of(x), sigma,
                                                    value_of(l)))>::value));
  EXPECT_TRUE((std::is_same<
               Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>,
               decltype(stan::math::gp_matern52_cov(
                   value_of(x), value_of(sigma), l))>::value));
}
