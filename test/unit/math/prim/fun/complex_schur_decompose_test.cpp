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
      EXPECT_NEAR(real(x(i, j)), real(y(i, j)), 1e-6);
      EXPECT_NEAR(imag(x(i, j)), imag(y(i, j)), 1e-6);
    }
  }
}

TEST(primFun, complex_schur_decompose) {
  using c_t = std::complex<double>;
  using stan::math::complex_schur_decompose_t;
  using stan::math::complex_schur_decompose_u;

  // verify that A = U T U*

  Eigen::MatrixXd A(3, 3);
  A << 0, 2, 2, 0, 0, 2, 1, 0, 1;
  auto A_t = complex_schur_decompose_t(A);
  auto A_u = complex_schur_decompose_u(A);
  auto A_recovered = A_u * A_t * A_u.adjoint();
  expect_complex_mat_eq(stan::math::to_complex(A,0), A_recovered);


  Eigen::MatrixXcd B(3, 3);
  B << 0, 2, 2, 0, c_t(0,1), 2, 1, 0, 1;

  auto B_t = complex_schur_decompose_t(B);
  auto B_u = complex_schur_decompose_u(B);

  auto B_recovered = B_u * B_t * B_u.adjoint();

  expect_complex_mat_eq(B, B_recovered);
}
