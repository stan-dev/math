#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(ProbExponential, quantile_values) {
  using stan::math::exponential_qf;

  EXPECT_FLOAT_EQ(exponential_qf(0.1, 8), 0.0131700644572);
  EXPECT_FLOAT_EQ(exponential_qf(0.4, 7), 0.0729750891094);
  EXPECT_FLOAT_EQ(exponential_qf(0.1, 7), 0.0150515022368);

  Eigen::VectorXd p(2);
  p << 0.2, 0.6;

  Eigen::VectorXd beta(2);
  beta << 16, 2;

  Eigen::VectorXd res(2);
  res << 0.0318776501877, 0.1308986759820;

  EXPECT_MATRIX_FLOAT_EQ(exponential_qf(p, 7), res);

  res << 0.0222921839962, 0.1783374719694;
  EXPECT_MATRIX_FLOAT_EQ(exponential_qf(0.3, beta), res);

  res << 0.0139464719571, 0.4581453659371;
  EXPECT_MATRIX_FLOAT_EQ(exponential_qf(p, beta), res);
}
