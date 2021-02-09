#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, ub) {
  EXPECT_FLOAT_EQ(2.0 - exp(-1.0), stan::math::ub_constrain(-1.0, 2.0));
  EXPECT_THROW(
      stan::math::ub_constrain(1.7, std::numeric_limits<double>::infinity()),
      std::domain_error);
}

TEST(prob_transform, ub_vec) {
  Eigen::VectorXd input(2);
  input << -1.0, 1.1;
  Eigen::VectorXd ubv(2);
  ubv << 2.0, 3.0;
  double ub = 2.0;

  Eigen::VectorXd resv(2);
  resv << 2.0 - exp(-1.0), 3.0 - exp(1.1);
  Eigen::VectorXd res(2);
  res << 2.0 - exp(-1.0), 2.0 - exp(1.1);

  EXPECT_MATRIX_EQ(resv, stan::math::ub_constrain(input, ubv));
  EXPECT_MATRIX_EQ(res, stan::math::ub_constrain(input, ub));

  double lp = 0.0;
  EXPECT_MATRIX_EQ(resv, stan::math::ub_constrain(input, ubv, lp));
  EXPECT_EQ(input.sum(), lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(res, stan::math::ub_constrain(input, ub, lp));
  EXPECT_EQ(input.sum(), lp);
}

TEST(prob_transform, ub_j) {
  double lp = 15.0;
  EXPECT_FLOAT_EQ(2.0 - exp(-1.0), stan::math::ub_constrain(-1.0, 2.0, lp));
  EXPECT_FLOAT_EQ(15.0 - 1.0, lp);

  double lp2 = 1.87;
  EXPECT_THROW(stan::math::ub_constrain(
                   -5.2, std::numeric_limits<double>::infinity(), lp2),
               std::domain_error);
  EXPECT_FLOAT_EQ(1.87, lp2);
}
TEST(prob_transform, ub_f) {
  double y = 2.0;
  double U = 4.0;
  EXPECT_FLOAT_EQ(log(-(y - U)), stan::math::ub_free(2.0, 4.0));

  EXPECT_THROW(
      stan::math::ub_free(19.765, std::numeric_limits<double>::infinity()),
      std::domain_error);
}
TEST(prob_transform, ub_f_exception) {
  double ub = 4.0;
  EXPECT_THROW(stan::math::ub_free(ub + 0.01, ub), std::domain_error);
}
TEST(prob_transform, ub_rt) {
  double x = -1.0;
  double xc = stan::math::ub_constrain(x, 2.0);
  double xcf = stan::math::ub_free(xc, 2.0);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::ub_constrain(xcf, 2.0);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
