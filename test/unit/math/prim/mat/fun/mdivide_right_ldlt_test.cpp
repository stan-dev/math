#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, mdivide_right_ldlt_val) {
  stan::math::LDLT_factor<double, -1, -1> ldlt_Ad;
  stan::math::matrix_d Ad(2, 2);
  Eigen::Matrix<double, 1, Eigen::Dynamic> b(2);

  Ad << 2.0, 3.0, 3.0, 7.0;
  b << 5.0, 6.0;

  ldlt_Ad.compute(Ad);
  ASSERT_TRUE(ldlt_Ad.success());

  auto res = mdivide_right_ldlt(b, ldlt_Ad);

  EXPECT_EQ(b.cols(), res.cols());
  EXPECT_EQ(b.rows(), res.rows());

  EXPECT_NEAR(res(0), 3.4, 1e-12);
  EXPECT_NEAR(res(1), -0.6, 1e-12);
}
