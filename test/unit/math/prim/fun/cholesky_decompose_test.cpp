#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMatrixPrimMat, cholesky_decompose) {
  stan::math::matrix_d m0;
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;

  using stan::math::cholesky_decompose;

  EXPECT_NO_THROW(cholesky_decompose(m0));
  EXPECT_THROW_MSG(cholesky_decompose(m1), std::invalid_argument,
                   "Expecting a square matrix");
}

TEST(MathMatrixPrimMat, cholesky_decompose_exception) {
  stan::math::matrix_d m;

  m.resize(2, 2);
  m << 1.0, 2.0, 2.0, 3.0;
  EXPECT_THROW_MSG(stan::math::cholesky_decompose(m), std::domain_error,
                   "Matrix m is not positive definite");

  m.resize(0, 0);
  EXPECT_NO_THROW(stan::math::cholesky_decompose(m));

  m.resize(2, 3);
  EXPECT_THROW_MSG(stan::math::cholesky_decompose(m), std::invalid_argument,
                   "Expecting a square matrix");

  // not symmetric
  m.resize(2, 2);
  m << 1.0, 2.0, 3.0, 4.0;
  EXPECT_THROW_MSG(stan::math::cholesky_decompose(m), std::domain_error,
                   "is not symmetric");
}

TEST(MathMatrixPrimMat, cholesky_decompose_expressions) {
  // Test for https://github.com/stan-dev/math/issues/3198
  stan::math::matrix_d A(2, 3);
  A << 1, 2, 3, 4, 5, 6;

  stan::math::vector_d L_u(3);
  L_u << 1, 0, 0.5;

  auto L = stan::math::cholesky_corr_constrain(L_u, 3);

  EXPECT_NO_THROW(stan::math::cholesky_decompose(stan::math::multiply(
      stan::math::multiply(A, stan::math::multiply_lower_tri_self_transpose(L)),
      stan::math::transpose(A))));
}
