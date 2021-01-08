#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, lub) {
  EXPECT_FLOAT_EQ(2.0 + (5.0 - 2.0) * stan::math::inv_logit(-1.0),
                  stan::math::lub_constrain(-1.0, 2.0, 5.0));

  EXPECT_THROW(
      stan::math::lub_constrain(1.7, -std::numeric_limits<double>::infinity(),
                                +std::numeric_limits<double>::infinity()),
      std::domain_error);
  EXPECT_THROW(stan::math::lub_constrain(
                   1.8, 3.0, +std::numeric_limits<double>::infinity()),
               std::domain_error);
  EXPECT_THROW(stan::math::lub_constrain(
                   1.9, -std::numeric_limits<double>::infinity(), -12.5),
               std::domain_error);
}

TEST(prob_transform, lub_vec) {
  Eigen::VectorXd input(2);
  input << -1.0, 1.1;
  Eigen::VectorXd lbv(2);
  lbv << -1.0, 0.5;
  Eigen::VectorXd ubv(2);
  ubv << 2.0, 3.0;
  double lb = 1.0;
  double ub = 2.0;

  Eigen::VectorXd resvv(2);
  resvv << (2.0 + 1.0) * stan::math::inv_logit(-1.0) - 1.0,
      (3.0 - 0.5) * stan::math::inv_logit(1.1) + 0.5;
  Eigen::VectorXd ressv(2);
  ressv << (2.0 - 1.0) * stan::math::inv_logit(-1.0) + 1.0,
      (3.0 - 1.0) * stan::math::inv_logit(1.1) + 1.0;
  Eigen::VectorXd resvs(2);
  resvs << (2.0 + 1.0) * stan::math::inv_logit(-1.0) - 1.0,
      (2.0 - 0.5) * stan::math::inv_logit(1.1) + 0.5;
  Eigen::VectorXd res(2);
  res << (2.0 - 1.0) * stan::math::inv_logit(-1.0) + 1.0,
      (2.0 - 1.0) * stan::math::inv_logit(1.1) + 1.0;

  EXPECT_MATRIX_EQ(resvv, stan::math::lub_constrain(input, lbv, ubv));
  EXPECT_MATRIX_EQ(ressv, stan::math::lub_constrain(input, lb, ubv));
  EXPECT_MATRIX_EQ(resvs, stan::math::lub_constrain(input, lbv, ub));
  EXPECT_MATRIX_EQ(res, stan::math::lub_constrain(input, lb, ub));

  double lp = 0.0;
  EXPECT_MATRIX_EQ(resvv, stan::math::lub_constrain(input, lbv, ubv, lp));
  EXPECT_FLOAT_EQ(((ubv - lbv).array().log() + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(ressv, stan::math::lub_constrain(input, lb, ubv, lp));
  EXPECT_FLOAT_EQ(((ubv.array() - lb).log() + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(resvs, stan::math::lub_constrain(input, lbv, ub, lp));
  EXPECT_FLOAT_EQ(((ub - lbv.array()).log() + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(res, stan::math::lub_constrain(input, lb, ub, lp));
  EXPECT_FLOAT_EQ((std::log(ub - lb) + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
}

TEST(prob_transform, lub_j) {
  double lp = -17.0;
  double L = 2.0;
  double U = 5.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + (U - L) * stan::math::inv_logit(x),
                  stan::math::lub_constrain(x, L, U, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U - L) + log(stan::math::inv_logit(x))
                      + log(1.0 - stan::math::inv_logit(x)),
                  lp);

  double lp1 = -12.9;
  EXPECT_THROW(
      stan::math::lub_constrain(1.7, -std::numeric_limits<double>::infinity(),
                                +std::numeric_limits<double>::infinity(), lp1),
      std::domain_error);
  EXPECT_FLOAT_EQ(-12.9, lp1);

  double lp2 = -19.8;
  double lp2_expected = -19.8;
  EXPECT_THROW(stan::math::lub_constrain(
                   1.8, 3.0, +std::numeric_limits<double>::infinity(), lp2),
               std::domain_error);
  EXPECT_FLOAT_EQ(lp2_expected, lp2);

  double lp3 = -422;
  double lp3_expected = -422;
  EXPECT_THROW(stan::math::lub_constrain(
                   1.9, -std::numeric_limits<double>::infinity(), -12.5, lp3),
               std::domain_error);
  EXPECT_FLOAT_EQ(lp3_expected, lp3);
}
TEST(prob_transform, lubException) {
  using stan::math::lub_constrain;
  EXPECT_THROW(lub_constrain(5.0, 1.0, 1.0), std::domain_error);
  EXPECT_NO_THROW(lub_constrain(5.0, 1.0, 1.01));
  double lp = 12;
  EXPECT_THROW(lub_constrain(5.0, 1.0, 1.0, lp), std::domain_error);
  EXPECT_NO_THROW(lub_constrain(5.0, 1.0, 1.01, lp));
}
TEST(prob_transform, lub_f) {
  double L = -10.0;
  double U = 27.0;
  double y = 3.0;
  EXPECT_FLOAT_EQ(stan::math::logit((y - L) / (U - L)),
                  stan::math::lub_free(y, L, U));

  EXPECT_THROW(
      stan::math::lub_free(14.2, -std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity()),
      std::domain_error);
  EXPECT_THROW(stan::math::lub_free(
                   -18.3, -std::numeric_limits<double>::infinity(), 7.6),
               std::domain_error);
  EXPECT_THROW(stan::math::lub_free(763.9, -3122.2,
                                    std::numeric_limits<double>::infinity()),
               std::domain_error);
}
TEST(prob_transform, lub_f_exception) {
  double L = -10.0;
  double U = 27.0;
  EXPECT_THROW(stan::math::lub_free(L - 0.01, L, U), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(U + 0.01, L, U), std::domain_error);

  EXPECT_THROW(stan::math::lub_free((L + U) / 2, U, L), std::domain_error);
}
TEST(prob_transform, lub_rt) {
  double x = -1.0;
  double xc = stan::math::lub_constrain(x, 2.0, 4.0);
  double xcf = stan::math::lub_free(xc, 2.0, 4.0);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::lub_constrain(xcf, 2.0, 4.0);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
TEST(prob_transform, underflow) {
  EXPECT_EQ(0, stan::math::lub_constrain(-1000, 0, 1));
  double lp = 0;
  EXPECT_EQ(0, stan::math::lub_constrain(-1000, 0, 1, lp));
}
