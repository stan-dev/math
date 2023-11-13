#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(primFun, complex_schur_decompose_ut) {
  using c_t = std::complex<double>;
  using stan::math::complex_schur_decompose_t;
  using stan::math::complex_schur_decompose_u;

  // verify that A = U T U*

  Eigen::MatrixXd A(3, 3);
  A << 0, 2, 2, 0, 0, 2, 1, 0, 1;
  auto A_t = complex_schur_decompose_t(A);
  auto A_u = complex_schur_decompose_u(A);
  auto A_recovered = A_u * A_t * A_u.adjoint();
  EXPECT_MATRIX_NEAR(A, stan::math::get_real(A_recovered), 1e-8);

  Eigen::MatrixXcd B(3, 3);
  B << 0, 2, 2, 0, c_t(0, 1), 2, 1, 0, 1;

  auto B_t = complex_schur_decompose_t(B);
  auto B_u = complex_schur_decompose_u(B);

  auto B_recovered = B_u * B_t * B_u.adjoint();
  EXPECT_MATRIX_COMPLEX_NEAR(B, B_recovered, 1e-8);
}

TEST(primFun, complex_schur_decompose) {
  using c_t = std::complex<double>;
  using stan::math::complex_schur_decompose;

  // verify that A = U T U*

  Eigen::MatrixXd A(3, 3);
  A << 0, 2, 2, 0, 0, 2, 1, 0, 1;
  Eigen::MatrixXcd A_t, A_u;
  std::tie(A_u, A_t) = complex_schur_decompose(A);
  auto A_recovered = A_u * A_t * A_u.adjoint();
  EXPECT_MATRIX_NEAR(A, stan::math::get_real(A_recovered), 1e-8);

  Eigen::MatrixXcd B(3, 3);
  B << 0, 2, 2, 0, c_t(0, 1), 2, 1, 0, 1;

  Eigen::MatrixXcd B_t, B_u;
  std::tie(B_u, B_t) = complex_schur_decompose(B);

  auto B_recovered = B_u * B_t * B_u.adjoint();
  EXPECT_MATRIX_COMPLEX_NEAR(B, B_recovered, 1e-8);
}
