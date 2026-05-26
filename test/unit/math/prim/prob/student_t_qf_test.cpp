#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, student_t_qf_vals) {
  using stan::math::student_t_qf;

  Eigen::VectorXd p{{0.1, 0.2, 0.5, 0.8, 0.9}};

  Eigen::VectorXd res = student_t_qf(p, 4, 2, 3);

  EXPECT_FLOAT_EQ(student_t_qf(0.2, 4, 2, 3), -0.8228937317);
  EXPECT_FLOAT_EQ(student_t_qf(0.8, 4, 2, 3), 4.822893732);
  EXPECT_FLOAT_EQ(student_t_qf(0.8, 3, 0, 1), 0.9784723124);

  EXPECT_FLOAT_EQ(student_t_qf(0.5, 3, 0, 1), 0.0);
}

TEST(MathFunctions, student_t_qf_vec) {
  using stan::math::student_t_qf;

  Eigen::VectorXd p{{0.1, 0.2, 0.5, 0.8, 0.9}};

  Eigen::VectorXd nu{{4, 4, 3, 3, 3}};

  Eigen::VectorXd res(5);
  for (int i = 0; i < 5; ++i) {
    res[i] = student_t_qf(p[i], nu[i], 2, 3);
  }
  EXPECT_MATRIX_FLOAT_EQ(res, student_t_qf(p, nu, 2, 3));

  Eigen::RowVectorXd p_rowvec = p.transpose();
  Eigen::RowVectorXd nu_rowvec = nu.transpose();
  Eigen::RowVectorXd res_rowvec = res.transpose();

  EXPECT_MATRIX_FLOAT_EQ(res, student_t_qf(p_rowvec, nu_rowvec, 2, 3));
}

TEST(MathFunctions, student_t_qf_inf) {
  using stan::math::student_t_qf;
  long double log_p = std::numeric_limits<long double>::min();
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(student_t_qf(0, 2, 0, 1), -inf);
  EXPECT_EQ(student_t_qf(1, 2, 0, 1), inf);
}

TEST(MathFunctions, student_t_qf_nan) {
  using stan::math::student_t_qf;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(student_t_qf(nan, 5, 0, 1), std::domain_error);
  EXPECT_THROW(student_t_qf(2.1, 5, 0, 1), std::domain_error);
}
