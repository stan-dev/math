#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, student_t_qf) {
  using stan::math::student_t_qf;
  EXPECT_FLOAT_EQ(0.0, student_t_qf(0.5, 2));

  EXPECT_NEAR(519.666270120263562, student_t_qf(0.95, 0.31), 1.5e-10);
  EXPECT_NEAR(-31830.988607907, student_t_qf(0.00001, 1), 1.5e-10);
}

TEST(MathFunctions, Equal) {
  using stan::math::student_t_qf;

  double p[] = {0.000000001, 0.0000001, 0.00001,   0.001,      0.05, 0.15, 0.25,
                0.35,        0.45,      0.55,      0.65,       0.75, 0.85, 0.95,
                0.999,       0.99999,   0.9999999, 0.999999999};

  double exact[]
      = {-234.02761040612,
         -73.98575804777,
         -23.33218270083,
         -7.17318221978,
         -2.13184678633,
         -1.18956685244,
         -0.74069708411,
         -0.41416326009,
         -0.13383036711,
          0.13383036711,   
          0.41416326009,
          0.74069708411,
          1.18956685244,
          2.13184678633,
          7.17318221978,
          23.33218270086,
          73.98575805751,
          234.02761206091};

  int numValues = sizeof(p) / sizeof(double);

  for (int i = 0; i < numValues; ++i) {
    EXPECT_NEAR(exact[i], student_t_qf(p[i], 4), 1.5e-11);
  }
}

TEST(MathFunctions, student_t_qf_inf) {
  using stan::math::student_t_qf;
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_EQ(student_t_qf(0, 2), -inf);
  EXPECT_EQ(student_t_qf(1.0, 2), inf);
}

TEST(MathFunctions, student_t_qf_nan) {
  using stan::math::student_t_qf;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(student_t_qf(0.3, nan), std::domain_error);
  EXPECT_THROW(student_t_qf(nan, 0.3), std::domain_error);
  EXPECT_THROW(student_t_qf(-2.0, 0.3), std::domain_error);
  EXPECT_THROW(student_t_qf(2.0, -0.3), std::domain_error);
}
