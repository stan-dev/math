#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, student_t_qf_vals) {
  using stan::math::student_t_qf;

  EXPECT_FLOAT_EQ(student_t_qf(0.2, 4, 2, 3), -0.8228937317);
  EXPECT_FLOAT_EQ(student_t_qf(0.8, 4, 2, 3), 4.822893732);
  EXPECT_FLOAT_EQ(student_t_qf(0.8, 3, 0, 1), 0.9784723124);

  EXPECT_FLOAT_EQ(student_t_qf(0.5, 3, 0, 1), 0.0);
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
