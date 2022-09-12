#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>

template <typename T, typename U>
void expect_complex_mat_eq(const T& x, const U& y, double tol = 1e-8) {
  EXPECT_EQ(x.rows(), y.rows());
  EXPECT_EQ(x.cols(), y.cols());
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      EXPECT_FLOAT_EQ(real(x(i, j)), real(y(i, j)));
      EXPECT_FLOAT_EQ(imag(x(i, j)), imag(y(i, j)));
    }
  }
}

TEST(primFun, complex_schur_decompose) {
  using c_t = std::complex<double>;
  using stan::math::complex_schur_decompose_t;
  using stan::math::complex_schur_decompose_u;

  // reference answers calculated using scipy.linalg.schur with output="complex"

  Eigen::MatrixXd A(3, 3);
  A << 0, 2, 2, 0, 1, 2, 1, 0, 1;

  auto A_t = complex_schur_decompose_t(A);

  Eigen::Matrix3cd A_t_expected;
  A_t_expected << c_t(2.65896708, 0.), c_t(-1.22839825, 1.32378589),
      c_t(0.42590094, 1.51937377), c_t(0., 0.), c_t(-0.32948354, 0.80225456),
      c_t(-0.59877805, 0.56192148), c_t(0., 0.), c_t(0., 0.),
      c_t(-0.32948354, -0.80225456);

  expect_complex_mat_eq(A_t_expected, A_t);
}
