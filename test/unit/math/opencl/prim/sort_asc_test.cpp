#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(OpenCLPrim, sort_asc_error) {
  stan::math::matrix_cl<double> a(2, 3);
  stan::math::matrix_cl<int> b(2, 3);

  EXPECT_THROW(stan::math::sort_asc(a), std::invalid_argument);
  EXPECT_THROW(stan::math::sort_asc(b), std::invalid_argument);
}

TEST(OpenCLPrim, sort_asc_double_simple) {
  Eigen::VectorXd a(9);
  a << 1, 2, 7, 4, 5, 6, 8, 3, 9;

  stan::math::matrix_cl<double> a_cl(a);

  Eigen::VectorXd correct = stan::math::sort_asc(a);
  Eigen::VectorXd res = stan::math::from_matrix_cl(stan::math::sort_asc(a_cl));

  EXPECT_MATRIX_EQ(res, correct);
}

TEST(OpenCLPrim, sort_asc_double_large) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(37985);

  stan::math::matrix_cl<double> a_cl(a);

  Eigen::VectorXd correct = stan::math::sort_asc(a);
  Eigen::VectorXd res = stan::math::from_matrix_cl(stan::math::sort_asc(a_cl));

  EXPECT_MATRIX_EQ(res, correct);
}

TEST(OpenCLPrim, sort_asc_int_simple) {
  Eigen::VectorXi a(9);
  a << 1, 2, 7, 4, 5, 6, 8, 3, 9;

  stan::math::matrix_cl<int> a_cl(a);

  Eigen::VectorXi correct = stan::math::sort_asc(a);
  Eigen::VectorXi res = stan::math::from_matrix_cl(stan::math::sort_asc(a_cl));

  EXPECT_MATRIX_EQ(res, correct);
}

TEST(OpenCLPrim, sort_asc_int_large) {
  Eigen::VectorXi a = Eigen::VectorXi::Random(37985);

  stan::math::matrix_cl<int> a_cl(a);

  Eigen::VectorXi correct = stan::math::sort_asc(a);
  Eigen::VectorXi res = stan::math::from_matrix_cl(stan::math::sort_asc(a_cl));

  EXPECT_MATRIX_EQ(res, correct);
}

#endif
