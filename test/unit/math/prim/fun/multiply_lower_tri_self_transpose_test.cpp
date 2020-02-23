#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

using stan::math::matrix_d;

matrix_d generate_large_L_tri_mat() {
  Eigen::Matrix<double, -1, -1> x(1000, 1000);

  x(0) = 0.1;
  for (int i = 1; i < 10000; ++i)
    x(i) = x(i - 1) + 0.1123456 * 1e10;
  return x;
}

void test_multiply_lower_tri_self_transpose(const matrix_d& x) {
  using stan::math::multiply_lower_tri_self_transpose;
  matrix_d y = multiply_lower_tri_self_transpose(x);
  matrix_d xp = x;
  for (int m = 0; m < xp.rows(); ++m)
    for (int n = m + 1; n < xp.cols(); ++n)
      xp(m, n) = 0;

  matrix_d xxt = xp * xp.transpose();
  EXPECT_EQ(y.rows(), xxt.rows());
  EXPECT_EQ(y.cols(), xxt.cols());
  for (int m = 0; m < y.rows(); ++m)
    for (int n = 0; n < y.cols(); ++n)
      EXPECT_FLOAT_EQ(xxt(m, n), y(m, n));
}

TEST(MathMatrixPrimMat, multiply_lower_tri_self_transpose) {
  using stan::math::check_symmetric;
  using stan::math::multiply_lower_tri_self_transpose;
  static const char* function
      = "stan::math::multiply_lower_tri_self_transpose(%1%)";
  Eigen::Matrix<double, -1, -1> x1;
  test_multiply_lower_tri_self_transpose(x1);

  Eigen::Matrix<double, -1, -1> x2(1, 1);
  x2 << 3.0;
  test_multiply_lower_tri_self_transpose(x2);

  Eigen::Matrix<double, -1, -1> x3(2, 2);
  x3 << 1.0, 0.0, 2.0, 3.0;
  test_multiply_lower_tri_self_transpose(x3);

  Eigen::Matrix<double, -1, -1> x4(3, 3);
  x4 << 1.0, 0.0, 0.0, 2.0, 3.0, 0.0, 4.0, 5.0, 6.0;
  test_multiply_lower_tri_self_transpose(x4);

  Eigen::Matrix<double, -1, -1> x5(3, 3);
  x5 << 1.0, 0.0, 100.0, 2.0, 3.0, 0.0, 4.0, 5.0, 6.0;
  test_multiply_lower_tri_self_transpose(x5);

  Eigen::Matrix<double, -1, -1> x6(3, 2);
  x6 << 1.0, 0.0, 2.0, 3.0, 4.0, 5.0;
  test_multiply_lower_tri_self_transpose(x6);

  Eigen::Matrix<double, -1, -1> x7(2, 3);
  x7 << 1.0, 0.0, 0.0, 2.0, 3.0, 0.0;
  test_multiply_lower_tri_self_transpose(x7);

  Eigen::Matrix<double, -1, -1> x8 = generate_large_L_tri_mat();
  EXPECT_NO_THROW(check_symmetric(function, "Symmetric matrix",
                                  multiply_lower_tri_self_transpose(x8)));
}
