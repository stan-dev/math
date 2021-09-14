#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, positive) {
  EXPECT_FLOAT_EQ(exp(-1.0), stan::math::positive_constrain(-1.0));
}
TEST(prob_transform, positive_j) {
  double lp = 15.0;
  EXPECT_FLOAT_EQ(exp(-1.0), stan::math::positive_constrain(-1.0, lp));
  EXPECT_FLOAT_EQ(15.0 - 1.0, lp);
}
TEST(prob_transform, positive_f) {
  EXPECT_FLOAT_EQ(log(0.5), stan::math::positive_free(0.5));
}
TEST(prob_transform, positive_f_exception) {
  EXPECT_THROW(stan::math::positive_free(-1.0), std::domain_error);
}
TEST(prob_transform, positive_rt) {
  double x = -1.0;
  double xc = stan::math::positive_constrain(x);
  double xcf = stan::math::positive_free(xc);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::positive_constrain(xcf);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
TEST(prob_transform, positive_vectorized) {
  double lp = 0;
  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  std::vector<Eigen::VectorXd> x_vec = {x, x, x};
  std::vector<Eigen::VectorXd> y_vec
      = stan::math::positive_constrain<false>(x_vec, lp);
  std::vector<Eigen::VectorXd> x_free_vec = stan::math::positive_free(y_vec);
  for (int i = 0; i < x_vec.size(); ++i) {
    EXPECT_MATRIX_FLOAT_EQ(x_vec[i], x_free_vec[i]);
  }
}
