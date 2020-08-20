#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsMultiNormalPrec, NotVectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.inverse();
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_prec_log(y, mu, L));
}

TEST(ProbDistributionsMultiNormalPrec, Vectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  vector<Matrix<double, Dynamic, 1> > vec_y(2);
  vector<Matrix<double, 1, Dynamic> > vec_y_t(2);
  Matrix<double, Dynamic, 1> y(3);
  Matrix<double, 1, Dynamic> y_t(3);
  y << 2.0, -2.0, 11.0;
  vec_y[0] = y;
  vec_y_t[0] = y;
  y << 4.0, -2.0, 1.0;
  vec_y[1] = y;
  vec_y_t[1] = y;
  y_t = y;

  vector<Matrix<double, Dynamic, 1> > vec_mu(2);
  vector<Matrix<double, 1, Dynamic> > vec_mu_t(2);
  Matrix<double, Dynamic, 1> mu(3);
  Matrix<double, 1, Dynamic> mu_t(3);
  mu << 1.0, -1.0, 3.0;
  vec_mu[0] = mu;
  vec_mu_t[0] = mu;
  mu << 2.0, -1.0, 4.0;
  vec_mu[1] = mu;
  vec_mu_t[1] = mu;
  mu_t = mu;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 10.0, -3.0, 0.0, -3.0, 5.0, 0.0, 0.0, 0.0, 5.0;
  Sigma = Sigma.inverse();

  // y and mu vectorized
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_prec_log(vec_y, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_prec_log(vec_y_t, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_prec_log(vec_y, vec_mu_t, Sigma));
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_prec_log(vec_y_t, vec_mu_t, Sigma));

  // y vectorized
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_prec_log(vec_y, mu, Sigma));
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_prec_log(vec_y_t, mu, Sigma));
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_prec_log(vec_y, mu_t, Sigma));
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_prec_log(vec_y_t, mu_t, Sigma));

  // mu vectorized
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_prec_log(y, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_prec_log(y_t, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_prec_log(y, vec_mu_t, Sigma));
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_prec_log(y_t, vec_mu_t, Sigma));
}
/*
TEST(ProbDistributionsMultiNormalPrec, MultiNormalOneRow) {
  Matrix<double, Dynamic, Dynamic> y(1, 3);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.inverse();
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_prec_log(y, mu, L));
}

TEST(ProbDistributionsMultiNormalPrec, MultiNormalMultiRow) {
  Matrix<double, Dynamic, Dynamic> y(2, 3);
  y << 2.0, -2.0, 11.0,
       4.0, -4.0, 22.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.inverse();
  EXPECT_FLOAT_EQ(-54.2152, stan::math::multi_normal_prec_log(y, mu, L));
}
*/

TEST(ProbDistributionsMultiNormalPrec, WrongSize) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  vector<Matrix<double, Dynamic, 1> > y_3_3(3);
  vector<Matrix<double, Dynamic, 1> > y_3_1(3);
  vector<Matrix<double, Dynamic, 1> > y_3_2(3);
  vector<Matrix<double, Dynamic, 1> > y_1_3(1);
  vector<Matrix<double, Dynamic, 1> > y_2_3(2);
  Matrix<double, Dynamic, 1> y_3(3);
  Matrix<double, Dynamic, 1> y_2(2);
  Matrix<double, Dynamic, 1> y_1(1);
  y_3 << 2.0, -2.0, 11.0;
  y_2 << 2.0, -2.0;
  y_1 << 2.0;
  y_3_3[0] = y_3;
  y_3_3[1] = y_3;
  y_3_3[2] = y_3;
  y_3_1[0] = y_1;
  y_3_1[1] = y_1;
  y_3_1[2] = y_1;
  y_3_2[0] = y_2;
  y_3_2[1] = y_2;
  y_3_2[2] = y_2;
  y_1_3[0] = y_3;
  y_2_3[0] = y_3;
  y_2_3[1] = y_3;

  vector<Matrix<double, Dynamic, 1> > mu_3_3(3);
  Matrix<double, Dynamic, 1> mu_3(3);
  mu_3 << 2.0, -2.0, 11.0;
  mu_3_3[0] = mu_3;
  mu_3_3[1] = mu_3;
  mu_3_3[2] = mu_3;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 10.0, -3.0, 0.0, -3.0, 5.0, 0.0, 0.0, 0.0, 5.0;
  Sigma = Sigma.inverse();

  EXPECT_NO_THROW(stan::math::multi_normal_prec_lpdf(y_3_3, mu_3_3, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_prec_lpdf(y_3, mu_3_3, Sigma));

  EXPECT_THROW(stan::math::multi_normal_prec_lpdf(y_1_3, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_prec_lpdf(y_2_3, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_prec_lpdf(y_3_1, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_prec_lpdf(y_3_2, mu_3_3, Sigma),
               std::invalid_argument);
}
